use crate::pack::{
    TrustedPackSink, TrustedPackedRecord, pack_trusted_fastq, pack_trusted_fastq_read_direct_sink,
    pack_trusted_fastq_read_sink,
};
use crate::{FastaReader, FastqConfig, FastqError, FastqReader, Result};
use std::io::BufRead;

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

    pub fn observe_sequence_record(&mut self, name: &[u8], seq: &[u8]) {
        self.records += 1;
        self.bases += seq.len() as u64;
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

pub fn consume_fasta<R: std::io::Read>(reader: &mut FastaReader<R>) -> Result<StreamStats> {
    let mut stats = StreamStats::default();
    reader.visit_records(|record| {
        stats.observe_sequence_record(record.name(), record.seq());
        Ok(())
    })?;
    Ok(stats)
}

pub fn consume_fasta_stats<R: std::io::Read>(reader: R) -> Result<StreamStats> {
    let mut reader = std::io::BufReader::new(reader);
    let mut line = Vec::new();
    let mut stats = StreamStats::default();
    let mut current_name_len = None;
    let mut current_bases = 0_u64;
    let mut current_first = 0_u8;
    let mut saw_record = false;

    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        let trimmed = trim_fasta_line(&line);
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with(b">") {
            if let Some(name_len) = current_name_len.replace(trimmed.len() as u64) {
                observe_fasta_stream_record(&mut stats, name_len, current_bases, current_first);
                current_bases = 0;
                current_first = 0;
            }
            if trimmed.len() == 1
                || trimmed[1..]
                    .first()
                    .is_some_and(|byte| byte.is_ascii_whitespace())
            {
                return Err(FastqError::Format("empty FASTA id".into()));
            }
            saw_record = true;
            continue;
        }
        if current_name_len.is_none() {
            return Err(FastqError::Format(
                "FASTA record header must start with `>`".into(),
            ));
        }
        if current_bases == 0 {
            current_first = trimmed[0];
        }
        current_bases += trimmed.len() as u64;
    }

    if let Some(name_len) = current_name_len {
        observe_fasta_stream_record(&mut stats, name_len, current_bases, current_first);
    } else if saw_record {
        return Err(FastqError::Format("empty FASTA id".into()));
    }
    Ok(stats)
}

fn observe_fasta_stream_record(stats: &mut StreamStats, name_len: u64, bases: u64, first_base: u8) {
    stats.records += 1;
    stats.bases += bases;
    stats.name_bytes += name_len;
    stats.checksum = stats
        .checksum
        .wrapping_add(first_base as u64)
        .wrapping_mul(1_099_511_628_211)
        .wrapping_add(bases);
}

fn trim_fasta_line(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    line.strip_suffix(b"\r").unwrap_or(line)
}

pub fn consume_trusted_fastq_with_pack(input: &[u8]) -> Result<StreamStats> {
    let mut stats = StreamStats::default();
    pack_trusted_fastq(input, |record| {
        observe_trusted_packed_record(&mut stats, record);
        Ok(())
    })?;
    Ok(stats)
}

pub fn consume_trusted_fastq_read_with_pack<R: std::io::Read>(
    reader: R,
    config: FastqConfig,
) -> Result<StreamStats> {
    let mut sink = StreamStatsSink::default();
    pack_trusted_fastq_read_sink(reader, config, &mut sink)?;
    Ok(sink.stats)
}

pub fn consume_trusted_fastq_read_direct_with_pack<R: std::io::Read>(
    reader: R,
    config: FastqConfig,
) -> Result<StreamStats> {
    let mut sink = StreamStatsSink::default();
    pack_trusted_fastq_read_direct_sink(reader, config, &mut sink)?;
    Ok(sink.stats)
}

fn observe_trusted_packed_record(stats: &mut StreamStats, record: TrustedPackedRecord<'_>) {
    stats.observe_record(record.name, record.seq, record.qual);
    stats.checksum = stats
        .checksum
        .wrapping_add(record.summary.bases.canonical_bases() as u64)
        .wrapping_add(record.summary.qualities.sum_phred);
}

#[derive(Default)]
struct StreamStatsSink {
    stats: StreamStats,
}

impl TrustedPackSink for &mut StreamStatsSink {
    fn record(&mut self, record: TrustedPackedRecord<'_>) -> Result<()> {
        observe_trusted_packed_record(&mut self.stats, record);
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SyntheticPattern {
    Cyclic,
    Entropy,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SyntheticFastaLayout {
    TwoLine,
    Wrapped { width: usize },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SyntheticAlphabet {
    Dna,
    Protein,
}

pub fn synthetic_fastq(records: usize, read_len: usize) -> Vec<u8> {
    synthetic_fastq_with_pattern(records, read_len, SyntheticPattern::Cyclic)
}

pub fn synthetic_fastq_with_pattern(
    records: usize,
    read_len: usize,
    pattern: SyntheticPattern,
) -> Vec<u8> {
    let mut out = Vec::with_capacity(records.saturating_mul(read_len + 32));
    for i in 0..records {
        push_fastq_record(&mut out, b"r", i, None, read_len, 0, pattern);
    }
    out
}

pub fn synthetic_interleaved_fastq_with_pattern(
    pairs: usize,
    read_len: usize,
    pattern: SyntheticPattern,
) -> Vec<u8> {
    let mut out = Vec::with_capacity(pairs.saturating_mul((read_len + 36) * 2));
    for i in 0..pairs {
        push_fastq_record(&mut out, b"frag", i, Some(1), read_len, 0, pattern);
        push_fastq_record(&mut out, b"frag", i, Some(2), read_len, 1, pattern);
    }
    out
}

pub fn synthetic_paired_fastq_with_pattern(
    pairs: usize,
    read_len: usize,
    pattern: SyntheticPattern,
) -> (Vec<u8>, Vec<u8>) {
    let mut r1 = Vec::with_capacity(pairs.saturating_mul(read_len + 36));
    let mut r2 = Vec::with_capacity(pairs.saturating_mul(read_len + 36));
    for i in 0..pairs {
        push_fastq_record(&mut r1, b"frag", i, Some(1), read_len, 0, pattern);
        push_fastq_record(&mut r2, b"frag", i, Some(2), read_len, 1, pattern);
    }
    (r1, r2)
}

pub fn synthetic_fasta(records: usize, read_len: usize) -> Vec<u8> {
    synthetic_fasta_with_options(
        records,
        read_len,
        SyntheticPattern::Cyclic,
        SyntheticFastaLayout::TwoLine,
        SyntheticAlphabet::Dna,
    )
}

pub fn synthetic_fasta_with_options(
    records: usize,
    read_len: usize,
    pattern: SyntheticPattern,
    layout: SyntheticFastaLayout,
    alphabet: SyntheticAlphabet,
) -> Vec<u8> {
    let mut out = Vec::with_capacity(records.saturating_mul(read_len + 16));
    for i in 0..records {
        out.extend_from_slice(b">r");
        push_usize_decimal(i, &mut out);
        out.push(b'\n');
        match layout {
            SyntheticFastaLayout::TwoLine => {
                push_symbols(&mut out, i, 0, read_len, pattern, alphabet);
                out.push(b'\n');
            }
            SyntheticFastaLayout::Wrapped { width } => {
                push_wrapped_symbols(&mut out, i, 0, read_len, pattern, alphabet, width.max(1));
            }
        }
    }
    out
}

fn push_fastq_record(
    out: &mut Vec<u8>,
    prefix: &[u8],
    index: usize,
    mate: Option<u8>,
    read_len: usize,
    phase: usize,
    pattern: SyntheticPattern,
) {
    out.push(b'@');
    out.extend_from_slice(prefix);
    push_usize_decimal(index, out);
    if let Some(mate) = mate {
        out.push(b'/');
        out.push(b'0' + mate);
    }
    out.push(b'\n');

    push_symbols(out, index, phase, read_len, pattern, SyntheticAlphabet::Dna);
    out.extend_from_slice(b"\n+\n");
    push_qualities(out, index, phase, read_len, pattern);
    out.push(b'\n');
}

fn push_wrapped_symbols(
    out: &mut Vec<u8>,
    index: usize,
    phase: usize,
    read_len: usize,
    pattern: SyntheticPattern,
    alphabet: SyntheticAlphabet,
    width: usize,
) {
    let mut written = 0;
    while written < read_len {
        let chunk = (read_len - written).min(width);
        push_symbols(out, index + written, phase, chunk, pattern, alphabet);
        out.push(b'\n');
        written += chunk;
    }
}

fn push_symbols(
    out: &mut Vec<u8>,
    index: usize,
    phase: usize,
    read_len: usize,
    pattern: SyntheticPattern,
    alphabet: SyntheticAlphabet,
) {
    let symbols = match alphabet {
        SyntheticAlphabet::Dna => &b"ACGT"[..],
        SyntheticAlphabet::Protein => &b"ACDEFGHIKLMNPQRSTVWY"[..],
    };
    match pattern {
        SyntheticPattern::Cyclic => {
            for j in 0..read_len {
                out.push(symbols[(index + j + phase) % symbols.len()]);
            }
        }
        SyntheticPattern::Entropy => {
            let mut state = rng_seed(index, phase, 0xa076_1d64_78bd_642f);
            for _ in 0..read_len {
                out.push(symbols[(next_u64(&mut state) as usize) % symbols.len()]);
            }
        }
    }
}

fn push_qualities(
    out: &mut Vec<u8>,
    index: usize,
    phase: usize,
    read_len: usize,
    pattern: SyntheticPattern,
) {
    match pattern {
        SyntheticPattern::Cyclic => out.extend(std::iter::repeat_n(b'I', read_len)),
        SyntheticPattern::Entropy => {
            let mut state = rng_seed(index, phase, 0xe703_7ed1_a0b4_28db);
            for _ in 0..read_len {
                out.push(33 + (next_u64(&mut state) % 41) as u8);
            }
        }
    }
}

fn rng_seed(index: usize, phase: usize, salt: u64) -> u64 {
    salt ^ ((index as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15)) ^ ((phase as u64) << 32)
}

fn next_u64(state: &mut u64) -> u64 {
    *state ^= *state >> 12;
    *state ^= *state << 25;
    *state ^= *state >> 27;
    state.wrapping_mul(0x2545_f491_4f6c_dd1d)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pack::pack_bases_and_summarize_qualities_into;
    use crate::{FastaReader, FastqReader};

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
    fn synthetic_fasta_fixture_parses_to_expected_stats() {
        let input = synthetic_fasta(3, 4);
        let mut reader = FastaReader::new(&input[..]);
        let stats = consume_fasta(&mut reader).unwrap();
        assert_eq!(stats.records, 3);
        assert_eq!(stats.bases, 12);
        assert_eq!(stats.qualities, 0);
    }

    #[test]
    fn streaming_fasta_stats_match_reader_stream_stats() {
        let input = synthetic_fasta_with_options(
            7,
            23,
            SyntheticPattern::Entropy,
            SyntheticFastaLayout::Wrapped { width: 5 },
            SyntheticAlphabet::Dna,
        );
        let mut reader = FastaReader::new(&input[..]);
        let reference = consume_fasta(&mut reader).unwrap();
        let streamed = consume_fasta_stats(&input[..]).unwrap();

        assert_eq!(streamed, reference);
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
