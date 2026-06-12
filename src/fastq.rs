use std::io::Read;
use std::ops::Range;

use crate::error::{FastqError, Result};
use crate::scan::scan_newlines;

const DEFAULT_SLAB_SIZE: usize = 8 * 1024 * 1024;

#[derive(Debug, Clone)]
pub struct FastqConfig {
    pub slab_size: usize,
    pub validate: bool,
}

impl Default for FastqConfig {
    fn default() -> Self {
        Self {
            slab_size: DEFAULT_SLAB_SIZE,
            validate: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecordRef {
    pub name: Range<u32>,
    pub seq: Range<u32>,
    pub plus: Range<u32>,
    pub qual: Range<u32>,
}

#[derive(Debug, Clone, Copy)]
pub struct FastqRecord<'a> {
    bytes: &'a [u8],
    record: &'a RecordRef,
}

impl<'a> FastqRecord<'a> {
    pub fn name(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.name.clone())]
    }

    pub fn seq(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.seq.clone())]
    }

    pub fn plus(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.plus.clone())]
    }

    pub fn qual(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.qual.clone())]
    }
}

#[derive(Debug)]
pub struct FastqBatch<'a> {
    bytes: &'a [u8],
    records: &'a [RecordRef],
}

impl<'a> FastqBatch<'a> {
    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn bytes(&self) -> &'a [u8] {
        self.bytes
    }

    pub fn record_refs(&self) -> &'a [RecordRef] {
        self.records
    }

    pub fn records(&self) -> impl Iterator<Item = FastqRecord<'a>> + 'a {
        self.records.iter().map(|record| FastqRecord {
            bytes: self.bytes,
            record,
        })
    }
}

#[derive(Debug)]
pub struct FastqReader<R> {
    reader: R,
    config: FastqConfig,
    buf: Vec<u8>,
    len: usize,
    next_start: usize,
    eof: bool,
    newlines: Vec<usize>,
    records: Vec<RecordRef>,
}

impl<R: Read> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        Self::with_config(reader, FastqConfig::default())
    }

    pub fn with_config(reader: R, config: FastqConfig) -> Self {
        let slab_size = config.slab_size.max(1024);
        let buf = vec![0; slab_size];
        Self {
            reader,
            config: FastqConfig {
                slab_size,
                ..config
            },
            buf,
            len: 0,
            next_start: 0,
            eof: false,
            newlines: Vec::with_capacity(slab_size / 48),
            records: Vec::with_capacity(8192),
        }
    }

    pub fn next_batch(&mut self) -> Result<Option<FastqBatch<'_>>> {
        self.compact_carry();
        self.fill_slab()?;

        self.records.clear();
        scan_newlines(&self.buf[..self.len], &mut self.newlines);
        self.next_start = frame_records(
            &self.buf[..self.len],
            &self.newlines,
            self.eof,
            self.config.validate,
            &mut self.records,
        )?;

        if self.records.is_empty() {
            if self.len == 0 && self.eof {
                return Ok(None);
            }
            return Err(FastqError::RecordTooLarge {
                slab_size: self.config.slab_size,
            });
        }

        Ok(Some(FastqBatch {
            bytes: &self.buf[..self.len],
            records: &self.records,
        }))
    }

    fn compact_carry(&mut self) {
        if self.next_start == 0 {
            return;
        }
        if self.next_start >= self.len {
            self.len = 0;
            self.next_start = 0;
            return;
        }
        let carry = self.len - self.next_start;
        self.buf.copy_within(self.next_start..self.len, 0);
        self.len = carry;
        self.next_start = 0;
    }

    fn fill_slab(&mut self) -> Result<()> {
        while !self.eof && self.len < self.config.slab_size {
            let n = self
                .reader
                .read(&mut self.buf[self.len..self.config.slab_size])?;
            if n == 0 {
                self.eof = true;
                break;
            }
            self.len += n;
        }
        Ok(())
    }
}

fn frame_records(
    bytes: &[u8],
    newline_offsets: &[usize],
    eof: bool,
    validate: bool,
    records: &mut Vec<RecordRef>,
) -> Result<usize> {
    let has_final_line = eof
        && newline_offsets
            .last()
            .map_or(!bytes.is_empty(), |&nl| nl + 1 < bytes.len());
    let line_count = newline_offsets.len() + usize::from(has_final_line);
    if eof && !line_count.is_multiple_of(4) {
        return Err(FastqError::Format("truncated FASTQ record".into()));
    }
    let complete_lines = (line_count / 4) * 4;

    for i in (0..complete_lines).step_by(4) {
        let name = line_range(bytes, newline_offsets, i);
        let seq = line_range(bytes, newline_offsets, i + 1);
        let plus = line_range(bytes, newline_offsets, i + 2);
        let qual = line_range(bytes, newline_offsets, i + 3);

        if validate {
            validate_record(bytes, &name, &seq, &plus, &qual)?;
        }

        records.push(RecordRef {
            name: to_u32(name)?,
            seq: to_u32(seq)?,
            plus: to_u32(plus)?,
            qual: to_u32(qual)?,
        });
    }

    if complete_lines == line_count {
        Ok(bytes.len())
    } else {
        Ok(line_start(newline_offsets, complete_lines))
    }
}

fn line_range(bytes: &[u8], newline_offsets: &[usize], line: usize) -> Range<usize> {
    let start = line_start(newline_offsets, line);
    let end = if line < newline_offsets.len() {
        newline_offsets[line]
    } else {
        bytes.len()
    };
    start..trim_cr_end(bytes, start, end)
}

fn line_start(newline_offsets: &[usize], line: usize) -> usize {
    if line == 0 {
        0
    } else {
        newline_offsets[line - 1] + 1
    }
}

fn validate_record(
    bytes: &[u8],
    name: &Range<usize>,
    seq: &Range<usize>,
    plus: &Range<usize>,
    qual: &Range<usize>,
) -> Result<()> {
    if bytes.get(name.start) != Some(&b'@') {
        return Err(FastqError::Format("header must start with `@`".into()));
    }
    if name.end == name.start + 1 {
        return Err(FastqError::Format("empty FASTQ id".into()));
    }
    if bytes.get(plus.start) != Some(&b'+') {
        return Err(FastqError::Format("plus line must start with `+`".into()));
    }
    let seq_len = seq.end - seq.start;
    let qual_len = qual.end - qual.start;
    if seq_len != qual_len {
        return Err(FastqError::Format(format!(
            "quality length {qual_len} != sequence length {seq_len}"
        )));
    }
    Ok(())
}

fn trim_cr_end(bytes: &[u8], start: usize, end: usize) -> usize {
    if end > start && bytes[end - 1] == b'\r' {
        end - 1
    } else {
        end
    }
}

fn to_u32(range: Range<usize>) -> Result<Range<u32>> {
    let start = u32::try_from(range.start)
        .map_err(|_| FastqError::Format("record offset exceeds u32 range".into()))?;
    let end = u32::try_from(range.end)
        .map_err(|_| FastqError::Format("record offset exceeds u32 range".into()))?;
    Ok(start..end)
}

fn to_usize(range: Range<u32>) -> Range<usize> {
    range.start as usize..range.end as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    fn collect_records(input: &[u8], slab_size: usize) -> Result<Vec<(Vec<u8>, Vec<u8>)>> {
        let mut reader = FastqReader::with_config(
            input,
            FastqConfig {
                slab_size,
                validate: true,
            },
        );
        let mut out = Vec::new();
        while let Some(batch) = reader.next_batch()? {
            for rec in batch.records() {
                out.push((rec.name().to_vec(), rec.seq().to_vec()));
            }
        }
        Ok(out)
    }

    #[test]
    fn reads_single_batch() {
        let out = collect_records(b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nJJJJ\n", 1024).unwrap();
        assert_eq!(
            out,
            vec![
                (b"@r1".to_vec(), b"ACGT".to_vec()),
                (b"@r2".to_vec(), b"TGCA".to_vec())
            ]
        );
    }

    #[test]
    fn handles_crlf() {
        let out = collect_records(b"@r1\r\nACGT\r\n+\r\nIIII\r\n", 1024).unwrap();
        assert_eq!(out, vec![(b"@r1".to_vec(), b"ACGT".to_vec())]);
    }

    #[test]
    fn accepts_missing_final_newline() {
        let out = collect_records(b"@r1\nACGT\n+\nIIII", 1024).unwrap();
        assert_eq!(out, vec![(b"@r1".to_vec(), b"ACGT".to_vec())]);
    }

    #[test]
    fn carries_split_records() {
        let out = collect_records(b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nJJJJ\n", 18).unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(out[1].1, b"TGCA");
    }

    #[test]
    fn rejects_bad_plus_line() {
        let err = collect_records(b"@r1\nACGT\n-\nIIII\n", 1024).unwrap_err();
        assert!(err.to_string().contains("plus line"));
    }

    #[test]
    fn rejects_quality_length_mismatch() {
        let err = collect_records(b"@r1\nACGT\n+\nIII\n", 1024).unwrap_err();
        assert!(
            err.to_string()
                .contains("quality length 3 != sequence length 4")
        );
    }

    #[test]
    fn rejects_truncated_eof() {
        let err = collect_records(b"@r1\nACGT\n+", 1024).unwrap_err();
        assert!(err.to_string().contains("truncated FASTQ record"));
    }
}
