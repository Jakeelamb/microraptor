use crate::{FastqReader, Result};

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

    #[test]
    fn synthetic_fixture_parses_to_expected_stats() {
        let input = synthetic_fastq(3, 4);
        let mut reader = FastqReader::new(&input[..]);
        let stats = consume_fastq(&mut reader).unwrap();
        assert_eq!(stats.records, 3);
        assert_eq!(stats.bases, 12);
        assert_eq!(stats.qualities, 12);
    }
}
