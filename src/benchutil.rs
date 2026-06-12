use crate::pack::{
    TrustedPackSink, TrustedPackedRecord, pack_trusted_fastq, pack_trusted_fastq_read_direct_sink,
    pack_trusted_fastq_read_sink,
};
use crate::{FastqConfig, FastqReader, Result};

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
