use std::io::Read;
use std::ops::Range;

use crate::error::{FastqError, Result};
use crate::fastq_frame::{self, RecordValidation};
use crate::scan::scan_newlines;

const DEFAULT_SLAB_SIZE: usize = 256 * 1024;

/// Configuration for FASTQ batch readers.
///
/// The default is a validated, unpaired reader with a 256 KiB slab and full pair
/// identifier validation when pairing APIs are used.
#[derive(Debug, Clone)]
pub struct FastqConfig {
    /// Target slab size in bytes.
    ///
    /// Values below 1024 are raised to 1024. Larger slabs reduce carry
    /// frequency for long records but increase the reusable buffer size.
    pub slab_size: usize,
    /// Validate FASTQ structure while framing records.
    ///
    /// Validation checks the leading `@`, leading `+`, and sequence/quality
    /// length equality. Disable only for trusted inputs.
    pub validate: bool,
    /// Whether a single stream should be interpreted as unpaired or interleaved.
    pub pairing: PairingMode,
    /// Identifier validation policy for paired APIs.
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
    /// Treat a single FASTQ stream as adjacent interleaved read pairs.
    pub fn interleaved(mut self) -> Self {
        self.pairing = PairingMode::Interleaved;
        self
    }

    /// Set the paired-read identifier validation policy.
    pub fn pair_validation(mut self, pair_validation: PairValidation) -> Self {
        self.pair_validation = pair_validation;
        self
    }
}

/// Pairing interpretation for a single FASTQ stream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairingMode {
    /// Records are yielded independently.
    None,
    /// Adjacent records are expected to be ordered mates.
    Interleaved,
}

/// Identifier validation policy for paired-end data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairValidation {
    /// Compare normalized first-token identifiers after removing `/1` or `/2`.
    Full,
    /// Fast path for ordered `/1` and `/2` mate suffixes, with fallback to
    /// [`Full`](Self::Full) when that shape is not present.
    FastSlash,
    /// Skip identifier checks for trusted, already synchronized inputs.
    None,
}

/// Byte ranges for a four-line FASTQ record within a batch slab.
///
/// Ranges are local to [`FastqBatch::bytes`]. They are exposed for callers that
/// want to build their own zero-copy views without using [`FastqRecord`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecordRef {
    /// Header line, including the leading `@` and excluding the newline.
    pub name: Range<u32>,
    /// Sequence line, excluding the newline.
    pub seq: Range<u32>,
    /// Plus line, excluding the newline.
    pub plus: Range<u32>,
    /// Quality line, excluding the newline.
    pub qual: Range<u32>,
}

/// Borrowed view of one FASTQ record inside a [`FastqBatch`].
///
/// The record borrows from the batch's slab buffer and cannot outlive the
/// batch. Accessors return byte slices rather than UTF-8 strings because FASTQ
/// names and qualities are byte-oriented.
#[derive(Debug, Clone, Copy)]
pub struct FastqRecord<'a> {
    bytes: &'a [u8],
    record: &'a RecordRef,
}

impl<'a> FastqRecord<'a> {
    /// Return the header line including the leading `@`.
    pub fn name(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.name.clone())]
    }

    /// Return the header line without a leading `@`.
    pub fn name_without_at(self) -> &'a [u8] {
        let name = self.name();
        name.strip_prefix(b"@").unwrap_or(name)
    }

    /// Return the first whitespace-delimited identifier token.
    pub fn id_token(self) -> &'a [u8] {
        let name = self.name_without_at();
        let end = name
            .iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(name.len());
        &name[..end]
    }

    /// Return the identifier token with a trailing `/1` or `/2` removed.
    pub fn pair_normalized_id(self) -> &'a [u8] {
        strip_pair_suffix(self.id_token())
    }

    /// Return the sequence line.
    pub fn seq(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.seq.clone())]
    }

    /// Return the plus line.
    pub fn plus(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.plus.clone())]
    }

    /// Return the quality line.
    pub fn qual(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.qual.clone())]
    }
}

/// Borrowed FASTQ record passed to [`FastqReader::visit_records`].
///
/// This view is optimized for single-pass consumers that do not need the batch
/// side table. The slices point directly into the reader slab and are valid
/// only for the duration of the visitor callback.
#[derive(Debug, Clone, Copy)]
pub struct FastqVisitRecord<'a> {
    name: &'a [u8],
    seq: &'a [u8],
    plus: &'a [u8],
    qual: &'a [u8],
}

impl<'a> FastqVisitRecord<'a> {
    /// Return the header line including the leading `@`.
    pub fn name(self) -> &'a [u8] {
        self.name
    }

    /// Return the sequence line.
    pub fn seq(self) -> &'a [u8] {
        self.seq
    }

    /// Return the plus line.
    pub fn plus(self) -> &'a [u8] {
        self.plus
    }

    /// Return the quality line.
    pub fn qual(self) -> &'a [u8] {
        self.qual
    }
}

/// Borrowed view of an ordered read pair.
#[derive(Debug, Clone, Copy)]
pub struct FastqPair<'a> {
    first: FastqRecord<'a>,
    second: FastqRecord<'a>,
}

impl<'a> FastqPair<'a> {
    /// First mate.
    pub fn first(&self) -> FastqRecord<'a> {
        self.first
    }

    /// Second mate.
    pub fn second(&self) -> FastqRecord<'a> {
        self.second
    }

    /// Normalized pair identifier from the first mate.
    pub fn pair_id(&self) -> &'a [u8] {
        self.first.pair_normalized_id()
    }
}

/// Strip a terminal `/1` or `/2` pair suffix from an identifier token.
pub fn strip_pair_suffix(id: &[u8]) -> &[u8] {
    fastq_frame::strip_pair_suffix(id)
}

/// A batch of borrowed FASTQ records from one reader slab.
///
/// The batch is invalidated by the next mutable call on the reader that
/// produced it. Process records before requesting another batch.
#[derive(Debug)]
pub struct FastqBatch<'a> {
    bytes: &'a [u8],
    records: &'a [RecordRef],
    base_offset: u64,
    first_record_index: u64,
    pair_validation: PairValidation,
}

impl<'a> FastqBatch<'a> {
    /// Number of records in the batch.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Whether the batch has no records.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Raw slab bytes backing this batch.
    pub fn bytes(&self) -> &'a [u8] {
        self.bytes
    }

    /// Record ranges within [`bytes`](Self::bytes).
    pub fn record_refs(&self) -> &'a [RecordRef] {
        self.records
    }

    /// Absolute byte offset of `bytes()[0]` in the original stream.
    pub fn base_offset(&self) -> u64 {
        self.base_offset
    }

    /// Zero-based index of the first record in this batch.
    pub fn first_record_index(&self) -> u64 {
        self.first_record_index
    }

    /// Pair validation mode inherited from the producing reader.
    pub fn pair_validation(&self) -> PairValidation {
        self.pair_validation
    }

    /// Iterate borrowed record views.
    pub fn records(&self) -> impl Iterator<Item = FastqRecord<'a>> + 'a {
        let bytes = self.bytes;
        let records: &'a [RecordRef] = self.records;
        records
            .iter()
            .map(move |record| FastqRecord { bytes, record })
    }

    /// Validate and iterate adjacent interleaved read pairs.
    ///
    /// Returns an error if the batch has an odd record count or mate
    /// identifiers fail the configured [`PairValidation`] mode.
    pub fn interleaved_pairs(&'a self) -> Result<InterleavedPairs<'a>> {
        validate_even_pair_count(self)?;
        validate_interleaved_pair_ids(self, self.pair_validation)?;
        Ok(InterleavedPairs {
            batch: self,
            next: 0,
        })
    }

    /// Validate and zip this batch with a mate batch from a separate reader.
    ///
    /// Identifier checks use this batch's configured [`PairValidation`] mode.
    pub fn paired_with(&'a self, mate: &'a FastqBatch<'a>) -> Result<PairedRecords<'a>> {
        validate_paired_batches(self, mate, self.pair_validation)?;
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

/// Iterator over adjacent pairs in an interleaved [`FastqBatch`].
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

/// Iterator over paired records from two separate [`FastqBatch`] values.
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

/// Validate and zip two separate FASTQ batches by ordered mate identifiers.
pub fn paired_records<'a>(
    first: &'a FastqBatch<'a>,
    second: &'a FastqBatch<'a>,
) -> Result<PairedRecords<'a>> {
    first.paired_with(second)
}

/// A paired-end batch produced by [`PairedFastqReader`].
///
/// The batch contains the validated common prefix from the two underlying
/// readers. If one reader yields more records than the other in a slab, the
/// extra records are retained for the next call.
#[derive(Debug)]
pub struct PairedFastqBatch<'a> {
    first_bytes: &'a [u8],
    first_records: &'a [RecordRef],
    second_bytes: &'a [u8],
    second_records: &'a [RecordRef],
}

impl<'a> PairedFastqBatch<'a> {
    /// Number of read pairs in the batch.
    pub fn len(&self) -> usize {
        self.first_records.len()
    }

    /// Whether this batch contains no read pairs.
    pub fn is_empty(&self) -> bool {
        self.first_records.is_empty()
    }

    /// Iterate borrowed read-pair views.
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

/// Iterator over read pairs in a [`PairedFastqBatch`].
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

/// Stateful reader for ordered separate-file paired-end FASTQ streams.
///
/// This reader validates matching records as it advances both streams. It does
/// not reorder or synchronize mates that appear in different orders.
#[derive(Debug)]
pub struct PairedFastqReader<R1, R2> {
    first: FastqReader<R1>,
    second: FastqReader<R2>,
    pair_validation: PairValidation,
}

impl<R1: Read, R2: Read> PairedFastqReader<R1, R2> {
    /// Create a paired reader with default FASTQ configuration.
    pub fn new(first: R1, second: R2) -> Self {
        Self::with_config(first, second, FastqConfig::default())
    }

    /// Create a paired reader using the same configuration for both streams.
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

    /// Create a paired reader with separate per-stream configurations.
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

    /// Create a paired reader from two existing [`FastqReader`] values.
    pub fn from_fastq_readers(first: FastqReader<R1>, second: FastqReader<R2>) -> Self {
        Self {
            pair_validation: first.config.pair_validation,
            first,
            second,
        }
    }

    /// Override the paired identifier validation policy.
    pub fn with_pair_validation(mut self, pair_validation: PairValidation) -> Self {
        self.pair_validation = pair_validation;
        self
    }

    /// Read the next validated paired batch.
    ///
    /// Returns `Ok(None)` when both streams end together.
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

/// Slab-based FASTQ reader over any [`Read`] input.
///
/// The reader frames four-line FASTQ records from a reusable byte slab. Records
/// crossing slab boundaries are carried into the next read. Returned batches
/// borrow from the reader and must be consumed before calling [`next_batch`]
/// again.
///
/// [`next_batch`]: Self::next_batch
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
    /// Create a reader with [`FastqConfig::default`].
    pub fn new(reader: R) -> Self {
        Self::with_config(reader, FastqConfig::default())
    }

    /// Create a reader with explicit configuration.
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
            newlines: Vec::new(),
            records: Vec::new(),
        }
    }

    /// Return the wrapped reader.
    pub fn into_inner(self) -> R {
        self.reader
    }

    /// Read the next batch of FASTQ records.
    ///
    /// Returns `Ok(None)` at EOF. Format errors include byte offset, record
    /// index, and line index when the failing location is known.
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
            pair_validation: self.config.pair_validation,
        }))
    }

    /// Visit every record in the stream without building a batch side table.
    ///
    /// This is the fastest public parse-only path for single-pass consumers.
    /// It still honors [`FastqConfig::validate`] for record structure, but it
    /// yields individual records and does not perform paired-end identifier
    /// validation. Use [`next_batch`](Self::next_batch),
    /// [`interleaved_pairs`](FastqBatch::interleaved_pairs), or
    /// [`PairedFastqReader`] when pair validation is required.
    pub fn visit_records<F>(&mut self, mut visit: F) -> Result<()>
    where
        F: FnMut(FastqVisitRecord<'_>) -> Result<()>,
    {
        loop {
            self.compact_carry();
            self.fill_slab()?;

            #[cfg(feature = "simd")]
            let (next_start, records) = {
                scan_newlines(&self.buf[..self.len], &mut self.newlines);
                visit_records_in_slab_from_newlines(
                    &self.buf[..self.len],
                    &self.newlines,
                    self.eof,
                    self.config.validate,
                    self.base_offset,
                    self.record_index,
                    &mut visit,
                )?
            };
            #[cfg(not(feature = "simd"))]
            let (next_start, records) = visit_records_in_slab(
                &self.buf[..self.len],
                self.eof,
                self.config.validate,
                self.base_offset,
                self.record_index,
                &mut visit,
            )?;
            self.next_start = next_start;
            self.record_index += records;

            if records == 0 {
                if self.len == 0 && self.eof {
                    return Ok(());
                }
                return Err(FastqError::RecordTooLarge {
                    slab_size: self.config.slab_size,
                });
            }
        }
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
            return Err(fastq_frame::format_at(
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

/// Visit records from an already resident FASTQ byte slice.
///
/// This path is intended for memory-mapped files, cached datasets, and other
/// callers that already own a complete FASTQ byte buffer. It validates the same
/// record structure as [`FastqReader::visit_records`] when
/// [`FastqConfig::validate`] is enabled, but it does not copy the input into a
/// streaming slab and does not perform paired-end identifier validation.
///
/// Returns the number of visited records.
pub fn visit_fastq_bytes<F>(bytes: &[u8], config: FastqConfig, mut visit: F) -> Result<u64>
where
    F: FnMut(FastqVisitRecord<'_>) -> Result<()>,
{
    #[cfg(feature = "simd")]
    let records = {
        let mut newlines = Vec::with_capacity(bytes.len() / 48);
        scan_newlines(bytes, &mut newlines);
        let (_, records) = visit_records_in_slab_from_newlines(
            bytes,
            &newlines,
            true,
            config.validate,
            0,
            0,
            &mut visit,
        )?;
        records
    };

    #[cfg(not(feature = "simd"))]
    let records = {
        let (_, records) = visit_records_in_slab(bytes, true, config.validate, 0, 0, &mut visit)?;
        records
    };

    Ok(records)
}

fn validate_even_pair_count(batch: &FastqBatch<'_>) -> Result<()> {
    if batch.records.len().is_multiple_of(2) {
        return Ok(());
    }
    let Some(last) = batch.records.last() else {
        return Ok(());
    };
    Err(fastq_frame::format_at(
        "interleaved FASTQ batch has an odd record count",
        batch.base_offset,
        last.name.start as usize,
        batch.first_record_index + (batch.records.len() - 1) as u64,
        0,
    ))
}

fn validate_paired_batches(
    first: &FastqBatch<'_>,
    second: &FastqBatch<'_>,
    pair_validation: PairValidation,
) -> Result<()> {
    if first.records.len() == second.records.len() {
        return validate_pair_ids(first, second, pair_validation);
    }

    let (batch, index) = if first.records.len() > second.records.len() {
        (first, second.records.len())
    } else {
        (second, first.records.len())
    };
    let record = &batch.records[index];
    Err(fastq_frame::format_at(
        "paired FASTQ batches have different record counts",
        batch.base_offset,
        record.name.start as usize,
        batch.first_record_index + index as u64,
        0,
    ))
}

fn validate_pair_ids(
    first: &FastqBatch<'_>,
    second: &FastqBatch<'_>,
    pair_validation: PairValidation,
) -> Result<()> {
    if pair_validation == PairValidation::None {
        return Ok(());
    }
    for index in 0..first.records.len() {
        let r1 = first.record_at(index);
        let r2 = second.record_at(index);
        if !pair_ids_match(r1, r2, pair_validation) {
            return Err(pair_id_mismatch(second, index));
        }
    }
    Ok(())
}

fn validate_interleaved_pair_ids(
    batch: &FastqBatch<'_>,
    pair_validation: PairValidation,
) -> Result<()> {
    if pair_validation == PairValidation::None {
        return Ok(());
    }
    for index in (0..batch.records.len()).step_by(2) {
        let r1 = batch.record_at(index);
        let r2 = batch.record_at(index + 1);
        if !pair_ids_match(r1, r2, pair_validation) {
            return Err(pair_id_mismatch(batch, index + 1));
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
        PairValidation::FastSlash => {
            fastq_frame::fast_slash_pair_ids_match(first.name(), second.name())
                .unwrap_or_else(|| first.pair_normalized_id() == second.pair_normalized_id())
        }
    }
}

fn pair_id_mismatch(batch: &FastqBatch<'_>, index: usize) -> FastqError {
    let record = &batch.records[index];
    fastq_frame::format_at(
        "paired FASTQ record identifiers do not match",
        batch.base_offset,
        record.name.start as usize,
        batch.first_record_index + index as u64,
        0,
    )
}

fn extra_record_error(batch: &FastqBatch<'_>, index: usize) -> FastqError {
    let record = &batch.records[index];
    fastq_frame::format_at(
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
    let layout = fastq_frame::slab_line_layout(bytes, newline_offsets, eof);
    if eof && !layout.line_count.is_multiple_of(4) {
        let start = fastq_frame::line_start(newline_offsets, (layout.line_count / 4) * 4);
        return Err(fastq_frame::format_at(
            "truncated FASTQ record",
            base_offset,
            start,
            first_record_index + records.len() as u64,
            (layout.line_count % 4) as u8,
        ));
    }
    if bytes.len() > u32::MAX as usize {
        return Err(FastqError::Format(
            "FASTQ slab byte offsets exceed u32 range".into(),
        ));
    }
    records.reserve(layout.complete_lines / 4);

    for i in (0..layout.complete_lines).step_by(4) {
        let record = fastq_frame::record_lines(bytes, newline_offsets, i);

        if validate {
            fastq_frame::validate_record(
                record,
                base_offset,
                first_record_index + records.len() as u64,
                RecordValidation::DEFAULT,
            )?;
        }

        records.push(RecordRef {
            name: to_u32_range(record.name.range()),
            seq: to_u32_range(record.seq.range()),
            plus: to_u32_range(record.plus.range()),
            qual: to_u32_range(record.qual.range()),
        });
    }

    if layout.complete_lines == layout.line_count && !layout.has_partial_trailing_line {
        Ok(bytes.len())
    } else {
        Ok(fastq_frame::line_start(
            newline_offsets,
            layout.complete_lines,
        ))
    }
}

#[cfg(feature = "simd")]
fn visit_records_in_slab_from_newlines<F>(
    bytes: &[u8],
    newline_offsets: &[usize],
    eof: bool,
    validate: bool,
    base_offset: u64,
    first_record_index: u64,
    visit: &mut F,
) -> Result<(usize, u64)>
where
    F: FnMut(FastqVisitRecord<'_>) -> Result<()>,
{
    let layout = fastq_frame::slab_line_layout(bytes, newline_offsets, eof);
    if eof && !layout.line_count.is_multiple_of(4) {
        let records = (layout.line_count / 4) as u64;
        return Err(fastq_frame::format_at(
            "truncated FASTQ record",
            base_offset,
            fastq_frame::line_start(newline_offsets, (layout.line_count / 4) * 4),
            first_record_index + records,
            (layout.line_count % 4) as u8,
        ));
    }
    let mut records = 0_u64;

    for line in (0..layout.complete_lines).step_by(4) {
        let record = fastq_frame::record_lines(bytes, newline_offsets, line);
        let record_index = first_record_index + records;

        if validate {
            fastq_frame::validate_record(
                record,
                base_offset,
                record_index,
                RecordValidation::DEFAULT,
            )?;
        }

        visit(FastqVisitRecord {
            name: record.name.bytes,
            seq: record.seq.bytes,
            plus: record.plus.bytes,
            qual: record.qual.bytes,
        })?;
        records += 1;
    }

    let next_start =
        if layout.complete_lines == layout.line_count && !layout.has_partial_trailing_line {
            bytes.len()
        } else {
            fastq_frame::line_start(newline_offsets, layout.complete_lines)
        };
    Ok((next_start, records))
}

#[cfg(not(feature = "simd"))]
fn visit_records_in_slab<F>(
    bytes: &[u8],
    eof: bool,
    validate: bool,
    base_offset: u64,
    first_record_index: u64,
    visit: &mut F,
) -> Result<(usize, u64)>
where
    F: FnMut(FastqVisitRecord<'_>) -> Result<()>,
{
    let mut cursor = 0;
    let mut records = 0_u64;
    let mut newlines = memchr::memchr_iter(b'\n', bytes);

    while cursor < bytes.len() {
        let record_start = cursor;
        let Some(name) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, eof) else {
            return Ok((record_start, records));
        };
        let Some(seq) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, eof) else {
            return incomplete_or_truncated_visit(
                eof,
                base_offset,
                record_start,
                first_record_index + records,
                records,
                1,
            );
        };
        let Some(plus) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, eof) else {
            return incomplete_or_truncated_visit(
                eof,
                base_offset,
                record_start,
                first_record_index + records,
                records,
                2,
            );
        };
        let Some(qual) = next_visit_line_bounds(bytes, &mut cursor, &mut newlines, eof) else {
            return incomplete_or_truncated_visit(
                eof,
                base_offset,
                record_start,
                first_record_index + records,
                records,
                3,
            );
        };

        let record_index = first_record_index + records;
        if validate {
            validate_visit_record_bounds(bytes, name, seq, plus, qual, base_offset, record_index)?;
        }

        visit(FastqVisitRecord {
            name: &bytes[name.0..name.1],
            seq: &bytes[seq.0..seq.1],
            plus: &bytes[plus.0..plus.1],
            qual: &bytes[qual.0..qual.1],
        })?;
        records += 1;
    }

    Ok((bytes.len(), records))
}

#[cfg(not(feature = "simd"))]
fn next_visit_line_bounds(
    bytes: &[u8],
    cursor: &mut usize,
    newlines: &mut impl Iterator<Item = usize>,
    eof: bool,
) -> Option<(usize, usize)> {
    let start = *cursor;
    if let Some(end) = newlines.next() {
        *cursor = end + 1;
        let end = fastq_frame::trim_cr_end(bytes, start, end);
        return Some((start, end));
    }
    if eof && start < bytes.len() {
        *cursor = bytes.len();
        let end = fastq_frame::trim_cr_end(bytes, start, bytes.len());
        return Some((start, end));
    }
    None
}

#[cfg(not(feature = "simd"))]
fn validate_visit_record_bounds(
    bytes: &[u8],
    name: (usize, usize),
    seq: (usize, usize),
    plus: (usize, usize),
    qual: (usize, usize),
    base_offset: u64,
    record_index: u64,
) -> Result<()> {
    if bytes.get(name.0) != Some(&b'@') {
        return Err(fastq_frame::format_at(
            "header must start with `@`",
            base_offset,
            name.0,
            record_index,
            0,
        ));
    }
    let id = &bytes[name.0 + 1..name.1];
    let id_end = id
        .iter()
        .position(u8::is_ascii_whitespace)
        .unwrap_or(id.len());
    if id_end == 0 {
        return Err(fastq_frame::format_at(
            "empty FASTQ id",
            base_offset,
            name.0,
            record_index,
            0,
        ));
    }
    if bytes.get(plus.0) != Some(&b'+') {
        return Err(fastq_frame::format_at(
            "plus line must start with `+`",
            base_offset,
            plus.0,
            record_index,
            2,
        ));
    }
    let seq_len = seq.1 - seq.0;
    let qual_len = qual.1 - qual.0;
    if seq_len != qual_len {
        return Err(fastq_frame::format_at(
            format!("quality length {qual_len} != sequence length {seq_len}"),
            base_offset,
            qual.0,
            record_index,
            3,
        ));
    }
    Ok(())
}

#[cfg(not(feature = "simd"))]
fn incomplete_or_truncated_visit(
    eof: bool,
    base_offset: u64,
    record_start: usize,
    record_index: u64,
    records: u64,
    line_index: u8,
) -> Result<(usize, u64)> {
    if eof {
        return Err(fastq_frame::format_at(
            "truncated FASTQ record",
            base_offset,
            record_start,
            record_index,
            line_index,
        ));
    }
    Ok((record_start, records))
}

fn to_u32_range(range: Range<usize>) -> Range<u32> {
    debug_assert!(u32::try_from(range.end).is_ok());
    range.start as u32..range.end as u32
}

fn to_usize(range: Range<u32>) -> Range<usize> {
    range.start as usize..range.end as usize
}

#[cfg(test)]
mod tests;
