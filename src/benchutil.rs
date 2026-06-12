use crate::pack::pack_bases_and_summarize_qualities_into;
use crate::scan::scan_newlines;
use crate::{FastqError, FastqPosition, FastqReader, Result};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct StreamStats {
    pub records: u64,
    pub bases: u64,
    pub qualities: u64,
    pub name_bytes: u64,
    pub checksum: u64,
}

impl StreamStats {
    pub fn observe_record(&mut self, name: &[u8], seq: &[u8], qual: &[u8]) {
        self.records += 1;
        self.bases += seq.len() as u64;
        self.qualities += qual.len() as u64;
        self.name_bytes += name.len() as u64;
        self.checksum = self
            .checksum
            .wrapping_add(seq.first().copied().unwrap_or_default() as u64)
            .wrapping_mul(1_099_511_628_211)
            .wrapping_add(seq.len() as u64);
    }
}

pub fn consume_fastq<R: std::io::Read>(reader: &mut FastqReader<R>) -> Result<StreamStats> {
    let mut stats = StreamStats::default();
    while let Some(batch) = reader.next_batch()? {
        for record in batch.records() {
            stats.observe_record(record.name(), record.seq(), record.qual());
        }
    }
    Ok(stats)
}

pub fn consume_trusted_fastq_with_pack(input: &[u8]) -> Result<StreamStats> {
    let mut stats = StreamStats::default();
    let mut packed = Vec::new();
    let mut mask = Vec::new();
    let mut newlines = Vec::with_capacity(input.len() / 48);
    scan_newlines(input, &mut newlines);

    let has_final_line = newlines
        .last()
        .map_or(!input.is_empty(), |&nl| nl + 1 < input.len());
    let line_count = newlines.len() + usize::from(has_final_line);
    if !line_count.is_multiple_of(4) {
        let record_index = (line_count / 4) as u64;
        return Err(format_at(
            "truncated FASTQ record",
            line_start(&newlines, (line_count / 4) * 4),
            record_index,
            (line_count % 4) as u8,
        ));
    }

    for line in (0..line_count).step_by(4) {
        let record_index = (line / 4) as u64;
        let name = line_range(input, &newlines, line);
        let seq = line_range(input, &newlines, line + 1);
        let plus = line_range(input, &newlines, line + 2);
        let qual = line_range(input, &newlines, line + 3);

        if input.get(name.start) != Some(&b'@') {
            return Err(format_at(
                "header must start with `@`",
                name.start,
                record_index,
                0,
            ));
        }
        if input.get(plus.start) != Some(&b'+') {
            return Err(format_at(
                "plus line must start with `+`",
                plus.start,
                record_index,
                2,
            ));
        }
        if seq.len() != qual.len() {
            return Err(format_at(
                format!(
                    "quality length {} != sequence length {}",
                    qual.len(),
                    seq.len()
                ),
                qual.start,
                record_index,
                3,
            ));
        }

        let summary =
            pack_bases_and_summarize_qualities_into(seq.bytes, qual.bytes, &mut packed, &mut mask)
                .map_err(|err| format_at(err.to_string(), qual.start, record_index, 3))?;
        stats.observe_record(name.bytes, seq.bytes, qual.bytes);
        stats.checksum = stats
            .checksum
            .wrapping_add(summary.bases.canonical_bases() as u64)
            .wrapping_add(summary.qualities.sum_phred);
    }

    Ok(stats)
}

pub fn synthetic_fastq(records: usize, read_len: usize) -> Vec<u8> {
    let bases = b"ACGT";
    let mut out = Vec::with_capacity(records.saturating_mul(read_len + 32));
    for i in 0..records {
        out.extend_from_slice(b"@r");
        push_usize_decimal(i, &mut out);
        out.push(b'\n');
        for j in 0..read_len {
            out.push(bases[(i + j) & 3]);
        }
        out.extend_from_slice(b"\n+\n");
        out.extend(std::iter::repeat_n(b'I', read_len));
        out.push(b'\n');
    }
    out
}

fn push_usize_decimal(mut n: usize, out: &mut Vec<u8>) {
    if n == 0 {
        out.push(b'0');
        return;
    }
    let mut buf = [0_u8; 20];
    let mut i = buf.len();
    while n != 0 {
        i -= 1;
        buf[i] = b'0' + (n % 10) as u8;
        n /= 10;
    }
    out.extend_from_slice(&buf[i..]);
}

#[derive(Clone, Copy)]
struct Line<'a> {
    bytes: &'a [u8],
    start: usize,
}

impl Line<'_> {
    fn len(self) -> usize {
        self.bytes.len()
    }
}

fn line_range<'a>(input: &'a [u8], newline_offsets: &[usize], line: usize) -> Line<'a> {
    let start = line_start(newline_offsets, line);
    let end = if line < newline_offsets.len() {
        newline_offsets[line]
    } else {
        input.len()
    };
    Line {
        bytes: &input[start..end],
        start,
    }
}

fn line_start(newline_offsets: &[usize], line: usize) -> usize {
    if line == 0 {
        0
    } else {
        newline_offsets[line - 1] + 1
    }
}

fn format_at(
    message: impl Into<String>,
    byte_offset: usize,
    record_index: u64,
    line_index: u8,
) -> FastqError {
    FastqError::FormatAt {
        message: message.into(),
        position: FastqPosition::new(byte_offset as u64, record_index, line_index),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::FastqReader;
    use crate::pack::pack_bases_and_summarize_qualities_into;

    #[test]
    fn synthetic_fixture_parses_to_expected_stats() {
        let input = synthetic_fastq(3, 4);
        let mut reader = FastqReader::new(&input[..]);
        let stats = consume_fastq(&mut reader).unwrap();
        assert_eq!(stats.records, 3);
        assert_eq!(stats.bases, 12);
        assert_eq!(stats.qualities, 12);
    }

    #[test]
    fn trusted_pack_matches_reader_pack_stats() {
        let input = synthetic_fastq(5, 7);
        let trusted = consume_trusted_fastq_with_pack(&input).unwrap();
        let mut reader = FastqReader::new(&input[..]);
        let mut reference = StreamStats::default();
        let mut packed = Vec::new();
        let mut mask = Vec::new();

        while let Some(batch) = reader.next_batch().unwrap() {
            for record in batch.records() {
                let summary = pack_bases_and_summarize_qualities_into(
                    record.seq(),
                    record.qual(),
                    &mut packed,
                    &mut mask,
                )
                .unwrap();
                reference.observe_record(record.name(), record.seq(), record.qual());
                reference.checksum = reference
                    .checksum
                    .wrapping_add(summary.bases.canonical_bases() as u64)
                    .wrapping_add(summary.qualities.sum_phred);
            }
        }

        assert_eq!(trusted, reference);
    }

    #[test]
    fn trusted_pack_accepts_missing_final_newline() {
        let input = b"@r0\nACGT\n+\nIIII";
        let stats = consume_trusted_fastq_with_pack(input).unwrap();
        assert_eq!(stats.records, 1);
        assert_eq!(stats.bases, 4);
    }

    #[test]
    fn trusted_pack_rejects_truncated_record() {
        let err = consume_trusted_fastq_with_pack(b"@r0\nACGT\n+\n").unwrap_err();
        assert!(err.to_string().contains("truncated FASTQ record"));
    }
}
