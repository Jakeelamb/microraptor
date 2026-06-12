use std::io::Read;
use std::ops::Range;

use crate::error::{FastqError, FastqPosition, Result};
use crate::scan::scan_newlines;

const DEFAULT_SLAB_SIZE: usize = 8 * 1024 * 1024;

#[derive(Debug, Clone)]
pub struct FastqConfig {
    pub slab_size: usize,
    pub validate: bool,
    pub pairing: PairingMode,
    pub pair_validation: PairValidation,
}

impl Default for FastqConfig {
    fn default() -> Self {
        Self {
            slab_size: DEFAULT_SLAB_SIZE,
            validate: true,
            pairing: PairingMode::None,
            pair_validation: PairValidation::Full,
        }
    }
}

impl FastqConfig {
    pub fn interleaved(mut self) -> Self {
        self.pairing = PairingMode::Interleaved;
        self
    }

    pub fn pair_validation(mut self, pair_validation: PairValidation) -> Self {
        self.pair_validation = pair_validation;
        self
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairingMode {
    None,
    Interleaved,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairValidation {
    Full,
    FastSlash,
    None,
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

    pub fn name_without_at(self) -> &'a [u8] {
        let name = self.name();
        name.strip_prefix(b"@").unwrap_or(name)
    }

    pub fn id_token(self) -> &'a [u8] {
        let name = self.name_without_at();
        let end = name
            .iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(name.len());
        &name[..end]
    }

    pub fn pair_normalized_id(self) -> &'a [u8] {
        strip_pair_suffix(self.id_token())
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

#[derive(Debug, Clone, Copy)]
pub struct FastqPair<'a> {
    first: FastqRecord<'a>,
    second: FastqRecord<'a>,
}

impl<'a> FastqPair<'a> {
    pub fn first(&self) -> FastqRecord<'a> {
        self.first
    }

    pub fn second(&self) -> FastqRecord<'a> {
        self.second
    }

    pub fn pair_id(&self) -> &'a [u8] {
        self.first.pair_normalized_id()
    }
}

pub fn strip_pair_suffix(id: &[u8]) -> &[u8] {
    if id.len() >= 2 && (id.ends_with(b"/1") || id.ends_with(b"/2")) {
        &id[..id.len() - 2]
    } else {
        id
    }
}

#[derive(Debug)]
pub struct FastqBatch<'a> {
    bytes: &'a [u8],
    records: &'a [RecordRef],
    base_offset: u64,
    first_record_index: u64,
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

    pub fn base_offset(&self) -> u64 {
        self.base_offset
    }

    pub fn first_record_index(&self) -> u64 {
        self.first_record_index
    }

    pub fn records(&self) -> impl Iterator<Item = FastqRecord<'a>> + 'a {
        let bytes = self.bytes;
        let records: &'a [RecordRef] = self.records;
        records
            .iter()
            .map(move |record| FastqRecord { bytes, record })
    }

    pub fn interleaved_pairs(&'a self) -> Result<InterleavedPairs<'a>> {
        validate_even_pair_count(self)?;
        validate_pair_ids(self, self)?;
        Ok(InterleavedPairs {
            batch: self,
            next: 0,
        })
    }

    pub fn paired_with(&'a self, mate: &'a FastqBatch<'a>) -> Result<PairedRecords<'a>> {
        validate_paired_batches(self, mate)?;
        Ok(PairedRecords {
            first: self,
            second: mate,
            next: 0,
        })
    }

    fn record_at(&self, index: usize) -> FastqRecord<'a> {
        let records: &'a [RecordRef] = self.records;
        FastqRecord {
            bytes: self.bytes,
            record: &records[index],
        }
    }
}

#[derive(Debug)]
pub struct InterleavedPairs<'a> {
    batch: &'a FastqBatch<'a>,
    next: usize,
}

impl<'a> Iterator for InterleavedPairs<'a> {
    type Item = FastqPair<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.batch.records.len() {
            return None;
        }
        let pair = FastqPair {
            first: self.batch.record_at(self.next),
            second: self.batch.record_at(self.next + 1),
        };
        self.next += 2;
        Some(pair)
    }
}

#[derive(Debug)]
pub struct PairedRecords<'a> {
    first: &'a FastqBatch<'a>,
    second: &'a FastqBatch<'a>,
    next: usize,
}

impl<'a> Iterator for PairedRecords<'a> {
    type Item = FastqPair<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.first.records.len() {
            return None;
        }
        let pair = FastqPair {
            first: self.first.record_at(self.next),
            second: self.second.record_at(self.next),
        };
        self.next += 1;
        Some(pair)
    }
}

pub fn paired_records<'a>(
    first: &'a FastqBatch<'a>,
    second: &'a FastqBatch<'a>,
) -> Result<PairedRecords<'a>> {
    first.paired_with(second)
}

#[derive(Debug)]
pub struct PairedFastqBatch<'a> {
    first_bytes: &'a [u8],
    first_records: &'a [RecordRef],
    second_bytes: &'a [u8],
    second_records: &'a [RecordRef],
}

impl<'a> PairedFastqBatch<'a> {
    pub fn len(&self) -> usize {
        self.first_records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.first_records.is_empty()
    }

    pub fn pairs(&'a self) -> PairedFastqPairs<'a> {
        PairedFastqPairs {
            batch: self,
            next: 0,
        }
    }

    fn first_record_at(&self, index: usize) -> FastqRecord<'a> {
        FastqRecord {
            bytes: self.first_bytes,
            record: &self.first_records[index],
        }
    }

    fn second_record_at(&self, index: usize) -> FastqRecord<'a> {
        FastqRecord {
            bytes: self.second_bytes,
            record: &self.second_records[index],
        }
    }
}

#[derive(Debug)]
pub struct PairedFastqPairs<'a> {
    batch: &'a PairedFastqBatch<'a>,
    next: usize,
}

impl<'a> Iterator for PairedFastqPairs<'a> {
    type Item = FastqPair<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.batch.first_records.len() {
            return None;
        }
        let pair = FastqPair {
            first: self.batch.first_record_at(self.next),
            second: self.batch.second_record_at(self.next),
        };
        self.next += 1;
        Some(pair)
    }
}

#[derive(Debug)]
pub struct PairedFastqReader<R1, R2> {
    first: FastqReader<R1>,
    second: FastqReader<R2>,
    pair_validation: PairValidation,
}

impl<R1: Read, R2: Read> PairedFastqReader<R1, R2> {
    pub fn new(first: R1, second: R2) -> Self {
        Self::with_config(first, second, FastqConfig::default())
    }

    pub fn with_config(first: R1, second: R2, config: FastqConfig) -> Self {
        Self::with_configs(
            first,
            FastqConfig {
                pairing: PairingMode::None,
                ..config.clone()
            },
            second,
            FastqConfig {
                pairing: PairingMode::None,
                ..config
            },
        )
    }

    pub fn with_configs(
        first: R1,
        first_config: FastqConfig,
        second: R2,
        second_config: FastqConfig,
    ) -> Self {
        let pair_validation = first_config.pair_validation;
        Self {
            first: FastqReader::with_config(
                first,
                FastqConfig {
                    pairing: PairingMode::None,
                    ..first_config
                },
            ),
            second: FastqReader::with_config(
                second,
                FastqConfig {
                    pairing: PairingMode::None,
                    ..second_config
                },
            ),
            pair_validation,
        }
    }

    pub fn from_fastq_readers(first: FastqReader<R1>, second: FastqReader<R2>) -> Self {
        Self {
            pair_validation: first.config.pair_validation,
            first,
            second,
        }
    }

    pub fn with_pair_validation(mut self, pair_validation: PairValidation) -> Self {
        self.pair_validation = pair_validation;
        self
    }

    pub fn next_pair_batch(&mut self) -> Result<Option<PairedFastqBatch<'_>>> {
        let first = self.first.next_batch()?;
        let second = self.second.next_batch()?;

        match (first, second) {
            (None, None) => Ok(None),
            (Some(first), None) => Err(extra_record_error(&first, 0)),
            (None, Some(second)) => Err(extra_record_error(&second, 0)),
            (Some(first), Some(second)) => {
                let pair_count = first.len().min(second.len());
                validate_pair_ids_prefix(&first, &second, pair_count, self.pair_validation)?;

                let first_view = BatchView::from_batch(&first, pair_count);
                let second_view = BatchView::from_batch(&second, pair_count);
                let first_retain = retain_from(&first, pair_count);
                let second_retain = retain_from(&second, pair_count);

                let _ = first;
                let _ = second;

                if let Some((next_start, record_count)) = first_retain {
                    self.first
                        .retain_records_from_parts(next_start, record_count);
                }
                if let Some((next_start, record_count)) = second_retain {
                    self.second
                        .retain_records_from_parts(next_start, record_count);
                }

                Ok(Some(PairedFastqBatch {
                    first_bytes: first_view.bytes(),
                    first_records: first_view.records(),
                    second_bytes: second_view.bytes(),
                    second_records: second_view.records(),
                }))
            }
        }
    }
}

struct BatchView {
    bytes_ptr: *const u8,
    bytes_len: usize,
    records_ptr: *const RecordRef,
    records_len: usize,
}

impl BatchView {
    fn from_batch(batch: &FastqBatch<'_>, records_len: usize) -> Self {
        Self {
            bytes_ptr: batch.bytes.as_ptr(),
            bytes_len: batch.bytes.len(),
            records_ptr: batch.records.as_ptr(),
            records_len,
        }
    }

    fn bytes<'a>(&self) -> &'a [u8] {
        // SAFETY: PairedFastqReader mutates only cursor/index fields between
        // saving this view and returning it. The backing buffer is not compacted
        // or reallocated until the next mutable reader call, which the returned
        // batch borrow prevents.
        unsafe { std::slice::from_raw_parts(self.bytes_ptr, self.bytes_len) }
    }

    fn records<'a>(&self) -> &'a [RecordRef] {
        // SAFETY: See bytes(). The saved length is capped to the validated
        // paired prefix and the records Vec is not cleared until the next
        // mutable reader call.
        unsafe { std::slice::from_raw_parts(self.records_ptr, self.records_len) }
    }
}

fn retain_from(batch: &FastqBatch<'_>, index: usize) -> Option<(usize, usize)> {
    if index >= batch.records.len() {
        return None;
    }
    Some((
        batch.records[index].name.start as usize,
        batch.records.len() - index,
    ))
}

#[derive(Debug)]
pub struct FastqReader<R> {
    reader: R,
    config: FastqConfig,
    buf: Vec<u8>,
    len: usize,
    next_start: usize,
    base_offset: u64,
    record_index: u64,
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
            base_offset: 0,
            record_index: 0,
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
        let first_record_index = self.record_index;
        self.next_start = frame_records(
            &self.buf[..self.len],
            &self.newlines,
            self.eof,
            self.config.validate,
            self.base_offset,
            first_record_index,
            &mut self.records,
        )?;
        self.align_interleaved_batch(first_record_index)?;

        self.record_index += self.records.len() as u64;

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
            base_offset: self.base_offset,
            first_record_index,
        }))
    }

    fn compact_carry(&mut self) {
        if self.next_start == 0 {
            return;
        }
        if self.next_start >= self.len {
            self.base_offset += self.len as u64;
            self.len = 0;
            self.next_start = 0;
            return;
        }
        let carry = self.len - self.next_start;
        self.buf.copy_within(self.next_start..self.len, 0);
        self.base_offset += self.next_start as u64;
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

    fn align_interleaved_batch(&mut self, first_record_index: u64) -> Result<()> {
        if self.config.pairing != PairingMode::Interleaved || self.records.len().is_multiple_of(2) {
            return Ok(());
        }

        let Some(last) = self.records.last() else {
            return Ok(());
        };
        if self.eof {
            return Err(format_at(
                "interleaved FASTQ ended with an unpaired record",
                self.base_offset,
                last.name.start as usize,
                first_record_index + (self.records.len() - 1) as u64,
                0,
            ));
        }

        self.next_start = last.name.start as usize;
        self.records.pop();
        Ok(())
    }

    fn retain_records_from_parts(&mut self, next_start: usize, record_count: usize) {
        self.next_start = next_start;
        self.record_index -= record_count as u64;
    }
}

fn validate_even_pair_count(batch: &FastqBatch<'_>) -> Result<()> {
    if batch.records.len().is_multiple_of(2) {
        return Ok(());
    }
    let Some(last) = batch.records.last() else {
        return Ok(());
    };
    Err(format_at(
        "interleaved FASTQ batch has an odd record count",
        batch.base_offset,
        last.name.start as usize,
        batch.first_record_index + (batch.records.len() - 1) as u64,
        0,
    ))
}

fn validate_paired_batches(first: &FastqBatch<'_>, second: &FastqBatch<'_>) -> Result<()> {
    if first.records.len() == second.records.len() {
        return validate_pair_ids(first, second);
    }

    let (batch, index) = if first.records.len() > second.records.len() {
        (first, second.records.len())
    } else {
        (second, first.records.len())
    };
    let record = &batch.records[index];
    Err(format_at(
        "paired FASTQ batches have different record counts",
        batch.base_offset,
        record.name.start as usize,
        batch.first_record_index + index as u64,
        0,
    ))
}

fn validate_pair_ids(first: &FastqBatch<'_>, second: &FastqBatch<'_>) -> Result<()> {
    if std::ptr::eq(first, second) {
        for index in (0..first.records.len()).step_by(2) {
            let r1 = first.record_at(index);
            let r2 = first.record_at(index + 1);
            if !pair_ids_match(r1, r2, PairValidation::Full) {
                return Err(pair_id_mismatch(first, index + 1));
            }
        }
        return Ok(());
    }

    for index in 0..first.records.len() {
        let r1 = first.record_at(index);
        let r2 = second.record_at(index);
        if !pair_ids_match(r1, r2, PairValidation::Full) {
            return Err(pair_id_mismatch(second, index));
        }
    }
    Ok(())
}

fn validate_pair_ids_prefix(
    first: &FastqBatch<'_>,
    second: &FastqBatch<'_>,
    len: usize,
    pair_validation: PairValidation,
) -> Result<()> {
    if pair_validation == PairValidation::None {
        return Ok(());
    }
    for index in 0..len {
        let r1 = first.record_at(index);
        let r2 = second.record_at(index);
        if !pair_ids_match(r1, r2, pair_validation) {
            return Err(pair_id_mismatch(second, index));
        }
    }
    Ok(())
}

fn pair_ids_match(first: FastqRecord<'_>, second: FastqRecord<'_>, mode: PairValidation) -> bool {
    match mode {
        PairValidation::None => true,
        PairValidation::Full => first.pair_normalized_id() == second.pair_normalized_id(),
        PairValidation::FastSlash => fast_slash_pair_ids_match(first.name(), second.name())
            .unwrap_or_else(|| first.pair_normalized_id() == second.pair_normalized_id()),
    }
}

fn fast_slash_pair_ids_match(first_name: &[u8], second_name: &[u8]) -> Option<bool> {
    let first = first_name.strip_prefix(b"@").unwrap_or(first_name);
    let second = second_name.strip_prefix(b"@").unwrap_or(second_name);
    if first.len() >= 3
        && first.len() == second.len()
        && first.ends_with(b"/1")
        && second.ends_with(b"/2")
    {
        return Some(first[..first.len() - 2] == second[..second.len() - 2]);
    }

    let first_end = token_end(first);
    let second_end = token_end(second);
    let first = &first[..first_end];
    let second = &second[..second_end];
    if first.len() < 3 || second.len() < 3 {
        return None;
    }
    let first_suffix = &first[first.len() - 2..];
    let second_suffix = &second[second.len() - 2..];
    if first_suffix != b"/1" || second_suffix != b"/2" || first.len() != second.len() {
        return None;
    }
    Some(first[..first.len() - 2] == second[..second.len() - 2])
}

fn token_end(bytes: &[u8]) -> usize {
    let mut end = 0;
    while end < bytes.len() && !bytes[end].is_ascii_whitespace() {
        end += 1;
    }
    end
}

fn pair_id_mismatch(batch: &FastqBatch<'_>, index: usize) -> FastqError {
    let record = &batch.records[index];
    format_at(
        "paired FASTQ record identifiers do not match",
        batch.base_offset,
        record.name.start as usize,
        batch.first_record_index + index as u64,
        0,
    )
}

fn extra_record_error(batch: &FastqBatch<'_>, index: usize) -> FastqError {
    let record = &batch.records[index];
    format_at(
        "paired FASTQ inputs have different record counts",
        batch.base_offset,
        record.name.start as usize,
        batch.first_record_index + index as u64,
        0,
    )
}

fn frame_records(
    bytes: &[u8],
    newline_offsets: &[usize],
    eof: bool,
    validate: bool,
    base_offset: u64,
    first_record_index: u64,
    records: &mut Vec<RecordRef>,
) -> Result<usize> {
    let has_final_line = eof
        && newline_offsets
            .last()
            .map_or(!bytes.is_empty(), |&nl| nl + 1 < bytes.len());
    let line_count = newline_offsets.len() + usize::from(has_final_line);
    if eof && !line_count.is_multiple_of(4) {
        let start = line_start(newline_offsets, (line_count / 4) * 4);
        return Err(format_at(
            "truncated FASTQ record",
            base_offset,
            start,
            first_record_index + records.len() as u64,
            (line_count % 4) as u8,
        ));
    }
    let complete_lines = (line_count / 4) * 4;

    for i in (0..complete_lines).step_by(4) {
        let name = line_range(bytes, newline_offsets, i);
        let seq = line_range(bytes, newline_offsets, i + 1);
        let plus = line_range(bytes, newline_offsets, i + 2);
        let qual = line_range(bytes, newline_offsets, i + 3);

        if validate {
            validate_record(
                bytes,
                &name,
                &seq,
                &plus,
                &qual,
                base_offset,
                first_record_index + records.len() as u64,
            )?;
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
    base_offset: u64,
    record_index: u64,
) -> Result<()> {
    if bytes.get(name.start) != Some(&b'@') {
        return Err(format_at(
            "header must start with `@`",
            base_offset,
            name.start,
            record_index,
            0,
        ));
    }
    if name.end == name.start + 1 {
        return Err(format_at(
            "empty FASTQ id",
            base_offset,
            name.start,
            record_index,
            0,
        ));
    }
    if bytes.get(plus.start) != Some(&b'+') {
        return Err(format_at(
            "plus line must start with `+`",
            base_offset,
            plus.start,
            record_index,
            2,
        ));
    }
    let seq_len = seq.end - seq.start;
    let qual_len = qual.end - qual.start;
    if seq_len != qual_len {
        return Err(format_at(
            format!("quality length {qual_len} != sequence length {seq_len}"),
            base_offset,
            qual.start,
            record_index,
            3,
        ));
    }
    Ok(())
}

fn format_at(
    message: impl Into<String>,
    base_offset: u64,
    local_offset: usize,
    record_index: u64,
    line_index: u8,
) -> FastqError {
    FastqError::FormatAt {
        message: message.into(),
        position: FastqPosition::new(base_offset + local_offset as u64, record_index, line_index),
    }
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
mod tests;
