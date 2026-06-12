use crate::pack::pack_bases_and_summarize_qualities_into;
use crate::scan::scan_newlines;
use crate::{FastqConfig, FastqError, FastqPosition, FastqReader, Result};

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
    consume_trusted_fastq_with_pack_into(input, 0, 0, true, &mut packed, &mut mask, &mut stats)?;
    Ok(stats)
}

pub fn consume_trusted_fastq_read_with_pack<R: std::io::Read>(
    mut reader: R,
    config: FastqConfig,
) -> Result<StreamStats> {
    let slab_size = config.slab_size.max(1024);
    let mut buf = vec![0_u8; slab_size];
    let mut len = 0;
    let mut eof = false;
    let mut base_offset = 0_u64;
    let mut record_index = 0_u64;
    let mut stats = StreamStats::default();
    let mut packed = Vec::new();
    let mut mask = Vec::new();

    loop {
        while !eof && len < slab_size {
            let n = reader.read(&mut buf[len..slab_size])?;
            if n == 0 {
                eof = true;
                break;
            }
            len += n;
        }

        let next_start = consume_trusted_fastq_with_pack_into(
            &buf[..len],
            base_offset,
            record_index,
            eof,
            &mut packed,
            &mut mask,
            &mut stats,
        )?;
        record_index = stats.records;

        if next_start == len {
            base_offset += len as u64;
            len = 0;
        } else {
            let carry = len - next_start;
            if next_start == 0 && carry == slab_size && !eof {
                return Err(FastqError::RecordTooLarge { slab_size });
            }
            buf.copy_within(next_start..len, 0);
            base_offset += next_start as u64;
            len = carry;
        }

        if eof {
            if len == 0 {
                return Ok(stats);
            }
            return Err(FastqError::RecordTooLarge { slab_size });
        }
    }
}

fn consume_trusted_fastq_with_pack_into(
    input: &[u8],
    base_offset: u64,
    first_record_index: u64,
    eof: bool,
    packed: &mut Vec<u8>,
    mask: &mut Vec<u8>,
    stats: &mut StreamStats,
) -> Result<usize> {
    let mut newlines = Vec::with_capacity(input.len() / 48);
    scan_newlines(input, &mut newlines);

    let has_final_line = eof
        && newlines
            .last()
            .map_or(!input.is_empty(), |&nl| nl + 1 < input.len());
    let line_count = newlines.len() + usize::from(has_final_line);
    if eof && !line_count.is_multiple_of(4) {
        let record_index = first_record_index + (line_count / 4) as u64;
        return Err(format_at(
            "truncated FASTQ record",
            base_offset,
            line_start(&newlines, (line_count / 4) * 4),
            record_index,
            (line_count % 4) as u8,
        ));
    }
    let complete_lines = (line_count / 4) * 4;

    for line in (0..complete_lines).step_by(4) {
        let record_index = first_record_index + (line / 4) as u64;
        let name = line_range(input, &newlines, line);
        let seq = line_range(input, &newlines, line + 1);
        let plus = line_range(input, &newlines, line + 2);
        let qual = line_range(input, &newlines, line + 3);

        observe_trusted_packed_record(
            name,
            seq,
            plus,
            qual,
            base_offset,
            record_index,
            packed,
            mask,
            stats,
        )?;
    }

    if complete_lines == line_count && eof {
        Ok(input.len())
    } else {
        let next_start = line_start(&newlines, complete_lines);
        if complete_lines == line_count && next_start == input.len() {
            return Ok(input.len());
        }
        Ok(next_start)
    }
}

#[allow(clippy::too_many_arguments)]
fn observe_trusted_packed_record(
    name: Line<'_>,
    seq: Line<'_>,
    plus: Line<'_>,
    qual: Line<'_>,
    base_offset: u64,
    record_index: u64,
    packed: &mut Vec<u8>,
    mask: &mut Vec<u8>,
    stats: &mut StreamStats,
) -> Result<()> {
    if name.bytes.first() != Some(&b'@') {
        return Err(format_at(
            "header must start with `@`",
            base_offset,
            name.start,
            record_index,
            0,
        ));
    }
    if plus.bytes.first() != Some(&b'+') {
        return Err(format_at(
            "plus line must start with `+`",
            base_offset,
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
            base_offset,
            qual.start,
            record_index,
            3,
        ));
    }

    let summary = pack_bases_and_summarize_qualities_into(seq.bytes, qual.bytes, packed, mask)
        .map_err(|err| format_at(err.to_string(), base_offset, qual.start, record_index, 3))?;
    stats.observe_record(name.bytes, seq.bytes, qual.bytes);
    stats.checksum = stats
        .checksum
        .wrapping_add(summary.bases.canonical_bases() as u64)
        .wrapping_add(summary.qualities.sum_phred);
    Ok(())
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
    let end = trim_cr_end(input, start, end);
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
    base_offset: u64,
    byte_offset: usize,
    record_index: u64,
    line_index: u8,
) -> FastqError {
    FastqError::FormatAt {
        message: message.into(),
        position: FastqPosition::new(base_offset + byte_offset as u64, record_index, line_index),
    }
}

fn trim_cr_end(bytes: &[u8], start: usize, end: usize) -> usize {
    if end > start && bytes[end - 1] == b'\r' {
        end - 1
    } else {
        end
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
        let reference = reference_pack_stats(&input, FastqConfig::default());

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

    #[test]
    fn trusted_stream_pack_matches_reader_across_slab_carry() {
        let input = synthetic_fastq(11, 37);
        let config = FastqConfig {
            slab_size: 128,
            ..FastqConfig::default()
        };
        let trusted = consume_trusted_fastq_read_with_pack(&input[..], config.clone()).unwrap();
        let reference = reference_pack_stats(&input, config);
        assert_eq!(trusted, reference);
    }

    #[test]
    fn trusted_stream_pack_trims_crlf_like_reader() {
        let input = b"@r0\r\nACGT\r\n+\r\nIIII\r\n@r1\r\nNN\r\n+\r\n!!\r\n";
        let config = FastqConfig {
            slab_size: 32,
            ..FastqConfig::default()
        };
        let trusted = consume_trusted_fastq_read_with_pack(&input[..], config.clone()).unwrap();
        let reference = reference_pack_stats(input, config);
        assert_eq!(trusted, reference);
    }

    #[test]
    fn trusted_stream_pack_accepts_missing_final_newline() {
        let input = b"@r0\nACGT\n+\nIIII";
        let config = FastqConfig {
            slab_size: 8,
            ..FastqConfig::default()
        };
        let stats = consume_trusted_fastq_read_with_pack(&input[..], config).unwrap();
        assert_eq!(stats.records, 1);
        assert_eq!(stats.bases, 4);
    }

    fn reference_pack_stats(input: &[u8], config: FastqConfig) -> StreamStats {
        let mut reader = FastqReader::with_config(input, config);
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
        reference
    }
}
