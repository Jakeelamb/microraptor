use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::ops::Range;

use crate::error::{FastqError, FastqPosition, Result};
#[cfg(feature = "bgzf")]
use crate::{BgzfDecodedBlockReader, BgzfIndex, BgzfIndexEntry, BgzfSeekReader, BgzfVirtualOffset};
use memchr::memchr;

const DEFAULT_BATCH_RECORDS: usize = 1024;
const DEFAULT_READER_BUFFER_SIZE: usize = 64 * 1024;
const REFERENCE_BATCH_RECORDS: usize = 16;
const REFERENCE_READER_BUFFER_SIZE: usize = 256 * 1024;
const REFERENCE_EXPECTED_SEQ_LEN: usize = 1024 * 1024;
const TWO_LINE_STREAM_BUFFER_SIZE: usize = 64 * 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ByteRange {
    start: u32,
    end: u32,
}

impl ByteRange {
    fn to_usize(self) -> Range<usize> {
        self.start as usize..self.end as usize
    }
}

/// Configuration for FASTA readers.
#[derive(Debug, Clone)]
pub struct FastaConfig {
    /// Maximum number of records returned per batch.
    ///
    /// Values below 1 are raised to 1. Larger batches reduce caller overhead
    /// but keep more sequence bytes resident until the next batch.
    pub batch_records: usize,
    /// Input buffer size used by streaming readers.
    ///
    /// Values below 1024 are raised to 1024. Larger buffers can reduce I/O
    /// calls for compressed or high-latency streams; smaller buffers may improve
    /// cache behavior on simple in-memory inputs.
    pub buffer_size: usize,
    /// Expected sequence length used as a preallocation hint.
    ///
    /// This is only a hint for reusable scratch buffers. It does not reject
    /// longer records.
    pub expected_seq_len: usize,
}

impl Default for FastaConfig {
    fn default() -> Self {
        Self {
            batch_records: DEFAULT_BATCH_RECORDS,
            buffer_size: DEFAULT_READER_BUFFER_SIZE,
            expected_seq_len: 0,
        }
    }
}

impl FastaConfig {
    /// Return parser settings tuned for chromosome-scale reference FASTA records.
    pub fn reference() -> Self {
        Self {
            batch_records: REFERENCE_BATCH_RECORDS,
            buffer_size: REFERENCE_READER_BUFFER_SIZE,
            expected_seq_len: REFERENCE_EXPECTED_SEQ_LEN,
        }
    }

    /// Return a copy configured for strict two-line FASTA-oriented streaming.
    pub fn two_line(mut self) -> Self {
        self.buffer_size = TWO_LINE_STREAM_BUFFER_SIZE;
        self
    }
}

/// Byte ranges for one FASTA record within a batch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecordRef {
    /// Header line, including the leading `>` and excluding the newline.
    pub name: Range<u32>,
    /// Sequence bytes with multiline FASTA sequence lines concatenated.
    pub seq: Range<u32>,
}

/// Borrowed view of one FASTA record inside a [`FastaBatch`].
#[derive(Debug, Clone, Copy)]
pub struct FastaRecord<'a> {
    bytes: &'a [u8],
    record: &'a FastaRecordRef,
    id_token: ByteRange,
}

impl<'a> FastaRecord<'a> {
    /// Return the header line including the leading `>`.
    pub fn name(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.name.clone())]
    }

    /// Return the header line without a leading `>`.
    pub fn name_without_gt(self) -> &'a [u8] {
        let name = self.name();
        name.strip_prefix(b">").unwrap_or(name)
    }

    /// Return the first whitespace-delimited identifier token.
    pub fn id_token(self) -> &'a [u8] {
        &self.bytes[self.id_token.to_usize()]
    }

    /// Return the concatenated sequence bytes.
    pub fn seq(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.seq.clone())]
    }
}

/// Borrowed FASTA record passed to [`FastaReader::visit_records`].
#[derive(Debug, Clone, Copy)]
pub struct FastaVisitRecord<'a> {
    name: &'a [u8],
    seq: &'a [u8],
}

impl<'a> FastaVisitRecord<'a> {
    /// Return the header line including the leading `>`.
    pub fn name(self) -> &'a [u8] {
        self.name
    }

    /// Return the header line without a leading `>`.
    pub fn name_without_gt(self) -> &'a [u8] {
        self.name.strip_prefix(b">").unwrap_or(self.name)
    }

    /// Return the first whitespace-delimited identifier token.
    pub fn id_token(self) -> &'a [u8] {
        let name = self.name_without_gt();
        let end = name
            .iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(name.len());
        &name[..end]
    }

    /// Return the concatenated sequence bytes.
    pub fn seq(self) -> &'a [u8] {
        self.seq
    }
}

const FASTA_STATS_CHECKSUM_INIT: u64 = 0xcbf2_9ce4_8422_2325;

/// Aggregate statistics for sequence-only FASTA workloads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FastaStats {
    /// Number of records observed.
    pub records: u64,
    /// Number of sequence bases observed.
    pub bases: u64,
    /// Lightweight deterministic checksum over record shape and edge bases.
    pub checksum: u64,
}

/// Detected FASTA physical layout for resident inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FastaShape {
    /// Empty input or only blank lines.
    Empty,
    /// Every non-empty record has exactly one header line and one sequence line.
    TwoLine,
    /// At least one record is multiline, blank-separated, or otherwise requires
    /// the robust FASTA parser.
    Multiline,
}

/// One `.fai`-style FASTA index entry.
///
/// `offset`, `line_bases`, and `line_width` follow the SAMtools `.fai`
/// convention over the uncompressed FASTA byte stream. When built from BGZF
/// input, `virtual_offset` stores the BGZF virtual offset for `offset`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaIndexEntry {
    /// Reference name, using the first whitespace-delimited token after `>`.
    pub name: Vec<u8>,
    /// Number of bases in the reference sequence.
    pub len: u64,
    /// Uncompressed byte offset of the first sequence byte.
    pub offset: u64,
    /// Number of bases per full sequence line.
    pub line_bases: u64,
    /// Number of bytes per full sequence line, including line ending bytes.
    pub line_width: u64,
    /// BGZF virtual offset for `offset`, when the index was built from BGZF.
    #[cfg(feature = "bgzf")]
    pub virtual_offset: Option<BgzfVirtualOffset>,
}

impl FastaIndexEntry {
    /// Return the uncompressed FASTA byte offset for a zero-based sequence
    /// position.
    ///
    /// Coordinates are 0-based and do not include FASTA line separators.
    pub fn sequence_offset(&self, pos: u64) -> Result<u64> {
        if pos > self.len {
            return Err(FastqError::Format(
                "FASTA sequence position exceeds reference length".into(),
            ));
        }
        if pos == self.len {
            return self.sequence_end_offset();
        }
        if self.line_bases == 0 {
            return Err(FastqError::Format(
                "FASTA index entry has zero line_bases for non-empty sequence".into(),
            ));
        }
        let line = pos / self.line_bases;
        let in_line = pos % self.line_bases;
        let line_offset = line.checked_mul(self.line_width).ok_or_else(|| {
            FastqError::Format("FASTA index offset calculation overflowed".into())
        })?;
        self.offset
            .checked_add(line_offset)
            .and_then(|offset| offset.checked_add(in_line))
            .ok_or_else(|| FastqError::Format("FASTA index offset calculation overflowed".into()))
    }

    /// Return physical FASTA byte spans covering a zero-based half-open
    /// sequence range.
    ///
    /// Each returned span points only at sequence bytes and excludes physical
    /// newline bytes.
    pub fn sequence_spans(&self, range: Range<u64>) -> Result<Vec<Range<u64>>> {
        self.validate_range(range.clone())?;
        let mut spans = Vec::new();
        let mut pos = range.start;
        while pos < range.end {
            if self.line_bases == 0 {
                return Err(FastqError::Format(
                    "FASTA index entry has zero line_bases for non-empty range".into(),
                ));
            }
            let in_line = pos % self.line_bases;
            let take = (range.end - pos).min(self.line_bases - in_line);
            let start = self.sequence_offset(pos)?;
            let end = start.checked_add(take).ok_or_else(|| {
                FastqError::Format("FASTA index span calculation overflowed".into())
            })?;
            spans.push(start..end);
            pos += take;
        }
        Ok(spans)
    }

    fn sequence_end_offset(&self) -> Result<u64> {
        if self.len == 0 {
            return Ok(self.offset);
        }
        self.sequence_offset(self.len - 1).map(|offset| offset + 1)
    }

    fn validate_range(&self, range: Range<u64>) -> Result<()> {
        if range.start > range.end {
            return Err(FastqError::Format(
                "FASTA range start must be <= end".into(),
            ));
        }
        if range.end > self.len {
            return Err(FastqError::Format(
                "FASTA range end exceeds reference length".into(),
            ));
        }
        Ok(())
    }
}

/// Configuration for planning balanced indexed-FASTA reference partitions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FastaPartitionConfig {
    /// Target number of partitions to produce.
    ///
    /// Values below 1 are raised to 1. Empty references are skipped.
    pub target_partitions: usize,
    /// Number of neighboring bases each partition should include in its
    /// fetch range on both sides of the core range.
    ///
    /// For k-mer or minimizer ingest, callers typically use the maximum
    /// lookaround needed by their k/window shape.
    pub overlap_bases: u64,
}

impl FastaPartitionConfig {
    /// Create a partition planner configuration.
    pub fn new(target_partitions: usize, overlap_bases: u64) -> Self {
        Self {
            target_partitions,
            overlap_bases,
        }
    }
}

/// One planned reference partition over indexed FASTA sequence coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaPartition {
    /// Zero-based partition index in output order.
    pub partition_index: usize,
    /// Reference name from the FASTA index.
    pub name: Vec<u8>,
    /// Core half-open sequence range assigned to this partition.
    pub core: Range<u64>,
    /// Half-open sequence range to fetch, expanded by overlap and clamped to
    /// the reference bounds.
    pub fetch: Range<u64>,
}

impl FastaPartition {
    /// Return the number of non-overlap bases assigned to this partition.
    pub fn core_len(&self) -> u64 {
        self.core.end - self.core.start
    }

    /// Return the zero-based offset of the first fetched base relative to the
    /// first core base.
    pub fn core_offset_in_fetch(&self) -> u64 {
        self.core.start - self.fetch.start
    }
}

/// Plan balanced reference partitions from a FASTA index.
///
/// The planner splits long references into approximately balanced chunks.
/// `core` ranges are disjoint and cover every non-empty reference; `fetch`
/// ranges include the requested overlap, clamped to each reference. The result
/// may contain more partitions than `target_partitions` when the index contains
/// many short references.
pub fn plan_fasta_partitions(
    index: &FastaIndex,
    config: FastaPartitionConfig,
) -> Result<Vec<FastaPartition>> {
    let target_partitions = config.target_partitions.max(1);
    let total_bases = index.entries().iter().try_fold(0_u64, |acc, entry| {
        acc.checked_add(entry.len)
            .ok_or_else(|| FastqError::Format("FASTA partition total length overflowed".into()))
    })?;
    if total_bases == 0 {
        return Ok(Vec::new());
    }
    let target_partitions_u64 = u64::try_from(target_partitions)
        .map_err(|_| FastqError::Format("FASTA partition count exceeds u64 range".into()))?;
    let target_bases = total_bases.div_ceil(target_partitions_u64).max(1);
    let mut partitions = Vec::new();

    for entry in index.entries() {
        let mut start = 0_u64;
        while start < entry.len {
            let remaining = entry.len - start;
            let take = remaining.min(target_bases);
            let core = start..start + take;
            let fetch_start = core.start.saturating_sub(config.overlap_bases);
            let fetch_end = core.end.saturating_add(config.overlap_bases).min(entry.len);
            partitions.push(FastaPartition {
                partition_index: partitions.len(),
                name: entry.name.clone(),
                core: core.clone(),
                fetch: fetch_start..fetch_end,
            });
            start = core.end;
        }
    }

    Ok(partitions)
}

/// Owned reference sequence chunk from an indexed FASTA source.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaReferenceChunk {
    /// Reference name from the FASTA index.
    pub name: Vec<u8>,
    /// Zero-based sequence offset of `seq[0]` within the reference.
    pub global_offset: u64,
    /// Sequence bytes for this chunk.
    pub seq: Vec<u8>,
}

/// Borrowed reference sequence chunk from an indexed FASTA source.
#[derive(Debug, Clone, Copy)]
pub struct FastaReferenceChunkRef<'a> {
    /// Reference name from the FASTA index.
    pub name: &'a [u8],
    /// Zero-based sequence offset of `seq[0]` within the reference.
    pub global_offset: u64,
    /// Sequence bytes for this chunk.
    pub seq: &'a [u8],
}

/// Sink trait for borrowed indexed FASTA reference chunks.
pub trait FastaReferenceChunkSink {
    /// Consume one borrowed reference chunk.
    fn chunk(&mut self, chunk: FastaReferenceChunkRef<'_>) -> Result<()>;
}

impl<F> FastaReferenceChunkSink for F
where
    F: FnMut(FastaReferenceChunkRef<'_>) -> Result<()>,
{
    fn chunk(&mut self, chunk: FastaReferenceChunkRef<'_>) -> Result<()> {
        self(chunk)
    }
}

/// Iterator over owned chunks from an [`IndexedFastaReader`].
pub struct FastaReferenceChunks<'a, R> {
    reader: &'a mut IndexedFastaReader<R>,
    name: Vec<u8>,
    next_offset: u64,
    end: u64,
    chunk_bases: u64,
}

/// Iterator over owned chunks from a [`BgzfIndexedFastaReader`].
#[cfg(feature = "bgzf")]
pub struct BgzfFastaReferenceChunks<'a, R> {
    reader: &'a mut BgzfIndexedFastaReader<R>,
    name: Vec<u8>,
    next_offset: u64,
    end: u64,
    chunk_bases: u64,
}

impl<R: Read + Seek> Iterator for FastaReferenceChunks<'_, R> {
    type Item = Result<FastaReferenceChunk>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next_offset >= self.end {
            return None;
        }
        let start = self.next_offset;
        let end = start.saturating_add(self.chunk_bases).min(self.end);
        self.next_offset = end;
        Some(
            self.reader
                .fetch(&self.name, start..end)
                .map(|seq| FastaReferenceChunk {
                    name: self.name.clone(),
                    global_offset: start,
                    seq,
                }),
        )
    }
}

#[cfg(feature = "bgzf")]
impl<R: Read + Seek> Iterator for BgzfFastaReferenceChunks<'_, R> {
    type Item = Result<FastaReferenceChunk>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next_offset >= self.end {
            return None;
        }
        let start = self.next_offset;
        let end = start.saturating_add(self.chunk_bases).min(self.end);
        self.next_offset = end;
        Some(
            self.reader
                .fetch(&self.name, start..end)
                .map(|seq| FastaReferenceChunk {
                    name: self.name.clone(),
                    global_offset: start,
                    seq,
                }),
        )
    }
}

/// A `.fai`-style FASTA index.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FastaIndex {
    entries: Vec<FastaIndexEntry>,
    name_to_index: HashMap<Vec<u8>, usize>,
}

impl FastaIndex {
    /// Return all index entries in FASTA order.
    pub fn entries(&self) -> &[FastaIndexEntry] {
        &self.entries
    }

    /// Return whether the index contains no references.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Return the number of references in the index.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Find an entry by reference name bytes.
    pub fn get(&self, name: &[u8]) -> Option<&FastaIndexEntry> {
        self.name_to_index
            .get(name)
            .and_then(|&idx| self.entries.get(idx))
    }

    /// Render the standard five-column `.fai` representation.
    pub fn to_fai_string(&self) -> String {
        let mut out = String::new();
        for entry in &self.entries {
            out.push_str(&String::from_utf8_lossy(&entry.name));
            out.push('\t');
            out.push_str(&entry.len.to_string());
            out.push('\t');
            out.push_str(&entry.offset.to_string());
            out.push('\t');
            out.push_str(&entry.line_bases.to_string());
            out.push('\t');
            out.push_str(&entry.line_width.to_string());
            out.push('\n');
        }
        out
    }

    /// Parse a standard five-column `.fai` index from bytes.
    pub fn from_fai_bytes(bytes: &[u8]) -> Result<Self> {
        let mut entries = Vec::new();
        let mut seen = HashSet::new();
        for (line_idx, line) in bytes.split(|&b| b == b'\n').enumerate() {
            let line = trim_line(line);
            if line.is_empty() {
                continue;
            }
            let fields = line.split(|&b| b == b'\t').collect::<Vec<_>>();
            if fields.len() != 5 {
                return Err(FastqError::Format(format!(
                    "invalid .fai line {}: expected 5 tab-delimited fields",
                    line_idx + 1
                )));
            }
            if fields[0].is_empty() {
                return Err(FastqError::Format(format!(
                    "invalid .fai line {}: empty reference name",
                    line_idx + 1
                )));
            }
            if !seen.insert(fields[0].to_vec()) {
                return Err(FastqError::Format(format!(
                    "invalid .fai line {}: duplicate reference name",
                    line_idx + 1
                )));
            }
            let len = parse_fai_u64(fields[1], line_idx + 1, "length")?;
            let offset = parse_fai_u64(fields[2], line_idx + 1, "offset")?;
            let line_bases = parse_fai_u64(fields[3], line_idx + 1, "line_bases")?;
            let line_width = parse_fai_u64(fields[4], line_idx + 1, "line_width")?;
            if len > 0 && line_bases == 0 {
                return Err(FastqError::Format(format!(
                    "invalid .fai line {}: non-empty reference has zero line_bases",
                    line_idx + 1
                )));
            }
            if line_width < line_bases {
                return Err(FastqError::Format(format!(
                    "invalid .fai line {}: line_width is smaller than line_bases",
                    line_idx + 1
                )));
            }
            entries.push(FastaIndexEntry {
                name: fields[0].to_vec(),
                len,
                offset,
                line_bases,
                line_width,
                #[cfg(feature = "bgzf")]
                virtual_offset: None,
            });
        }
        Ok(Self::from_entries(entries))
    }

    /// Parse a standard five-column `.fai` index from UTF-8 text.
    pub fn from_fai_str(text: &str) -> Result<Self> {
        Self::from_fai_bytes(text.as_bytes())
    }

    /// Parse a standard five-column `.fai` index from a reader.
    pub fn from_fai_read<R: Read>(mut reader: R) -> Result<Self> {
        let mut bytes = Vec::new();
        reader.read_to_end(&mut bytes)?;
        Self::from_fai_bytes(&bytes)
    }

    fn from_entries(entries: Vec<FastaIndexEntry>) -> Self {
        let name_to_index = entries
            .iter()
            .enumerate()
            .map(|(idx, entry)| (entry.name.clone(), idx))
            .collect();
        Self {
            entries,
            name_to_index,
        }
    }
}

fn parse_fai_u64(value: &[u8], line: usize, field: &str) -> Result<u64> {
    let text = std::str::from_utf8(value).map_err(|_| {
        FastqError::Format(format!(
            "invalid .fai line {line}: {field} is not valid UTF-8"
        ))
    })?;
    text.parse::<u64>().map_err(|_| {
        FastqError::Format(format!(
            "invalid .fai line {line}: {field} is not an unsigned integer"
        ))
    })
}

fn indexed_physical_window(entry: &FastaIndexEntry, range: Range<u64>) -> Result<Range<u64>> {
    entry.validate_range(range.clone())?;
    let start = entry.sequence_offset(range.start)?;
    let end = if range.start == range.end {
        start
    } else {
        entry
            .sequence_offset(range.end - 1)?
            .checked_add(1)
            .ok_or_else(|| FastqError::Format("FASTA index span calculation overflowed".into()))?
    };
    if end < start {
        return Err(FastqError::Format(
            "FASTA index physical range is invalid".into(),
        ));
    }
    Ok(start..end)
}

fn copy_indexed_sequence_window<R: Read>(
    reader: &mut R,
    entry: &FastaIndexEntry,
    mut physical_offset: u64,
    mut len: u64,
    expected_bases: usize,
    out: &mut Vec<u8>,
    scratch: &mut Vec<u8>,
) -> Result<()> {
    if len == 0 {
        return Ok(());
    }
    if entry.line_width == 0 {
        return Err(FastqError::Format(
            "FASTA index entry has zero line_width for non-empty range".into(),
        ));
    }
    if scratch.is_empty() {
        scratch.resize(8192, 0);
    }
    while len > 0 {
        let take = usize::try_from(len.min(scratch.len() as u64))
            .map_err(|_| FastqError::Format("FASTA fetch span exceeds usize range".into()))?;
        reader.read_exact(&mut scratch[..take])?;
        for &byte in &scratch[..take] {
            let rel = physical_offset.checked_sub(entry.offset).ok_or_else(|| {
                FastqError::Format("FASTA index physical offset precedes sequence offset".into())
            })?;
            if rel % entry.line_width < entry.line_bases {
                out.push(byte);
                if out.len() > expected_bases {
                    return Err(FastqError::Format(
                        "FASTA fetch produced more bases than expected".into(),
                    ));
                }
            }
            physical_offset = physical_offset.checked_add(1).ok_or_else(|| {
                FastqError::Format("FASTA fetch physical offset overflowed".into())
            })?;
        }
        len -= take as u64;
    }
    if out.len() != expected_bases {
        return Err(FastqError::Format(
            "FASTA fetch produced fewer bases than expected".into(),
        ));
    }
    Ok(())
}

/// Seekable FASTA reader backed by a `.fai` index.
pub struct IndexedFastaReader<R> {
    inner: R,
    index: FastaIndex,
    scratch: Vec<u8>,
}

impl<R: Read + Seek> IndexedFastaReader<R> {
    /// Create a seekable FASTA reader from an input stream and index.
    pub fn new(inner: R, index: FastaIndex) -> Self {
        Self {
            inner,
            index,
            scratch: Vec::new(),
        }
    }

    /// Return the loaded FASTA index.
    pub fn index(&self) -> &FastaIndex {
        &self.index
    }

    /// Fetch a zero-based half-open sequence range into `out`.
    pub fn fetch_into(&mut self, name: &[u8], range: Range<u64>, out: &mut Vec<u8>) -> Result<()> {
        let entry = self.index.get(name).ok_or_else(|| {
            FastqError::Format(format!(
                "FASTA reference not found in index: {}",
                String::from_utf8_lossy(name)
            ))
        })?;
        let window = indexed_physical_window(entry, range.clone())?;
        out.clear();
        let expected_bases = usize::try_from(range.end - range.start).map_err(|_| {
            FastqError::Format("FASTA fetch range length exceeds usize range".into())
        })?;
        out.reserve(expected_bases);
        self.inner.seek(SeekFrom::Start(window.start))?;
        copy_indexed_sequence_window(
            &mut self.inner,
            entry,
            window.start,
            window.end - window.start,
            expected_bases,
            out,
            &mut self.scratch,
        )?;
        Ok(())
    }

    /// Fetch a zero-based half-open sequence range into an owned buffer.
    pub fn fetch(&mut self, name: &[u8], range: Range<u64>) -> Result<Vec<u8>> {
        let mut out = Vec::new();
        self.fetch_into(name, range, &mut out)?;
        Ok(out)
    }

    /// Stream owned sequence chunks for one reference range.
    ///
    /// `chunk_bases` values below 1 are raised to 1. Chunks are yielded as
    /// owned buffers so callers can move them to worker threads after each
    /// iterator step.
    pub fn reference_chunks(
        &mut self,
        name: &[u8],
        range: Range<u64>,
        chunk_bases: u64,
    ) -> Result<FastaReferenceChunks<'_, R>> {
        let entry = self.index.get(name).ok_or_else(|| {
            FastqError::Format(format!(
                "FASTA reference not found in index: {}",
                String::from_utf8_lossy(name)
            ))
        })?;
        entry.validate_range(range.clone())?;
        Ok(FastaReferenceChunks {
            reader: self,
            name: name.to_vec(),
            next_offset: range.start,
            end: range.end,
            chunk_bases: chunk_bases.max(1),
        })
    }

    /// Stream borrowed sequence chunks into a sink using caller-owned storage.
    ///
    /// `chunk_bases` values below 1 are raised to 1. `out` is cleared and
    /// reused for each chunk, and the borrowed chunk passed to `sink` is valid
    /// only until the next sink call or until this method returns.
    pub fn reference_chunks_into<S>(
        &mut self,
        name: &[u8],
        range: Range<u64>,
        chunk_bases: u64,
        out: &mut Vec<u8>,
        sink: &mut S,
    ) -> Result<()>
    where
        S: FastaReferenceChunkSink,
    {
        let entry = self.index.get(name).ok_or_else(|| {
            FastqError::Format(format!(
                "FASTA reference not found in index: {}",
                String::from_utf8_lossy(name)
            ))
        })?;
        entry.validate_range(range.clone())?;
        stream_reference_chunks_into(name, range, chunk_bases, out, sink, |name, range, out| {
            self.fetch_into(name, range, out)
        })
    }

    /// Fetch a planned partition into an owned reference chunk.
    pub fn fetch_partition(&mut self, partition: &FastaPartition) -> Result<FastaReferenceChunk> {
        self.fetch(&partition.name, partition.fetch.clone())
            .map(|seq| FastaReferenceChunk {
                name: partition.name.clone(),
                global_offset: partition.fetch.start,
                seq,
            })
    }

    /// Return the wrapped reader and index.
    pub fn into_inner(self) -> (R, FastaIndex) {
        (self.inner, self.index)
    }
}

/// Seekable BGZF-compressed FASTA reader backed by `.fai` and BGZF block
/// indexes.
#[cfg(feature = "bgzf")]
pub struct BgzfIndexedFastaReader<R> {
    inner: BgzfSeekReader<R>,
    fasta_index: FastaIndex,
    bgzf_index: BgzfIndex,
    scratch: Vec<u8>,
}

#[cfg(feature = "bgzf")]
impl<R: Read + Seek> BgzfIndexedFastaReader<R> {
    /// Create a BGZF FASTA random-access reader.
    pub fn new(inner: R, fasta_index: FastaIndex, bgzf_index: BgzfIndex) -> Self {
        Self {
            inner: BgzfSeekReader::new(inner),
            fasta_index,
            bgzf_index,
            scratch: Vec::new(),
        }
    }

    /// Return the loaded FASTA index.
    pub fn fasta_index(&self) -> &FastaIndex {
        &self.fasta_index
    }

    /// Return the loaded BGZF index.
    pub fn bgzf_index(&self) -> &BgzfIndex {
        &self.bgzf_index
    }

    /// Fetch a zero-based half-open sequence range into `out`.
    pub fn fetch_into(&mut self, name: &[u8], range: Range<u64>, out: &mut Vec<u8>) -> Result<()> {
        let entry = self.fasta_index.get(name).ok_or_else(|| {
            FastqError::Format(format!(
                "FASTA reference not found in index: {}",
                String::from_utf8_lossy(name)
            ))
        })?;
        let window = indexed_physical_window(entry, range.clone())?;
        out.clear();
        let expected_bases = usize::try_from(range.end - range.start).map_err(|_| {
            FastqError::Format("FASTA fetch range length exceeds usize range".into())
        })?;
        out.reserve(expected_bases);
        if expected_bases == 0 {
            return Ok(());
        }
        let virtual_offset = self
            .bgzf_index
            .virtual_offset_for_uncompressed_offset(window.start)?
            .ok_or_else(|| FastqError::Bgzf("BGZF span offset is not indexed".into()))?;
        self.inner.seek_virtual_offset(virtual_offset)?;
        copy_indexed_sequence_window(
            &mut self.inner,
            entry,
            window.start,
            window.end - window.start,
            expected_bases,
            out,
            &mut self.scratch,
        )?;
        Ok(())
    }

    /// Fetch a zero-based half-open sequence range into an owned buffer.
    pub fn fetch(&mut self, name: &[u8], range: Range<u64>) -> Result<Vec<u8>> {
        let mut out = Vec::new();
        self.fetch_into(name, range, &mut out)?;
        Ok(out)
    }

    /// Stream owned sequence chunks for one reference range.
    ///
    /// `chunk_bases` values below 1 are raised to 1. Chunks are yielded as
    /// owned buffers so callers can move them to worker threads after each
    /// iterator step.
    pub fn reference_chunks(
        &mut self,
        name: &[u8],
        range: Range<u64>,
        chunk_bases: u64,
    ) -> Result<BgzfFastaReferenceChunks<'_, R>> {
        let entry = self.fasta_index.get(name).ok_or_else(|| {
            FastqError::Format(format!(
                "FASTA reference not found in index: {}",
                String::from_utf8_lossy(name)
            ))
        })?;
        entry.validate_range(range.clone())?;
        Ok(BgzfFastaReferenceChunks {
            reader: self,
            name: name.to_vec(),
            next_offset: range.start,
            end: range.end,
            chunk_bases: chunk_bases.max(1),
        })
    }

    /// Stream borrowed sequence chunks into a sink using caller-owned storage.
    ///
    /// `chunk_bases` values below 1 are raised to 1. `out` is cleared and
    /// reused for each chunk, and the borrowed chunk passed to `sink` is valid
    /// only until the next sink call or until this method returns.
    pub fn reference_chunks_into<S>(
        &mut self,
        name: &[u8],
        range: Range<u64>,
        chunk_bases: u64,
        out: &mut Vec<u8>,
        sink: &mut S,
    ) -> Result<()>
    where
        S: FastaReferenceChunkSink,
    {
        let entry = self.fasta_index.get(name).ok_or_else(|| {
            FastqError::Format(format!(
                "FASTA reference not found in index: {}",
                String::from_utf8_lossy(name)
            ))
        })?;
        entry.validate_range(range.clone())?;
        stream_reference_chunks_into(name, range, chunk_bases, out, sink, |name, range, out| {
            self.fetch_into(name, range, out)
        })
    }

    /// Fetch a planned partition into an owned reference chunk.
    pub fn fetch_partition(&mut self, partition: &FastaPartition) -> Result<FastaReferenceChunk> {
        self.fetch(&partition.name, partition.fetch.clone())
            .map(|seq| FastaReferenceChunk {
                name: partition.name.clone(),
                global_offset: partition.fetch.start,
                seq,
            })
    }
}

fn stream_reference_chunks_into<S, F>(
    name: &[u8],
    range: Range<u64>,
    chunk_bases: u64,
    out: &mut Vec<u8>,
    sink: &mut S,
    mut fetch_into: F,
) -> Result<()>
where
    S: FastaReferenceChunkSink,
    F: FnMut(&[u8], Range<u64>, &mut Vec<u8>) -> Result<()>,
{
    let chunk_bases = chunk_bases.max(1);
    let mut next_offset = range.start;
    while next_offset < range.end {
        let start = next_offset;
        let end = start.saturating_add(chunk_bases).min(range.end);
        next_offset = end;
        fetch_into(name, start..end, out)?;
        sink.chunk(FastaReferenceChunkRef {
            name,
            global_offset: start,
            seq: out,
        })?;
    }
    Ok(())
}

impl Default for FastaStats {
    fn default() -> Self {
        Self {
            records: 0,
            bases: 0,
            checksum: FASTA_STATS_CHECKSUM_INIT,
        }
    }
}

impl FastaStats {
    /// Observe one sequence.
    pub fn observe_sequence(&mut self, seq: &[u8]) {
        self.observe_sequence_parts(
            seq.len() as u64,
            seq.first().copied().unwrap_or_default(),
            seq.last().copied().unwrap_or_default(),
        );
    }

    fn observe_sequence_parts(&mut self, len: u64, first: u8, last: u8) {
        self.records += 1;
        self.bases += len;
        self.checksum ^= len;
        self.checksum = self
            .checksum
            .rotate_left(5)
            .wrapping_mul(0x0000_0100_0000_01b3);
        self.checksum ^= first as u64;
        self.checksum ^= (last as u64) << 8;
    }
}

/// Sink trait for FASTA record visitors.
pub trait FastaRecordSink {
    /// Consume one borrowed FASTA record.
    fn record(&mut self, record: FastaVisitRecord<'_>) -> Result<()>;
}

impl<F> FastaRecordSink for F
where
    F: FnMut(FastaVisitRecord<'_>) -> Result<()>,
{
    fn record(&mut self, record: FastaVisitRecord<'_>) -> Result<()> {
        self(record)
    }
}

/// A batch of borrowed FASTA records.
#[derive(Debug)]
pub struct FastaBatch<'a> {
    bytes: &'a [u8],
    records: &'a [FastaRecordRef],
    id_tokens: &'a [ByteRange],
    first_record_index: u64,
}

impl<'a> FastaBatch<'a> {
    /// Number of records in the batch.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Whether the batch has no records.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Raw bytes backing this batch.
    pub fn bytes(&self) -> &'a [u8] {
        self.bytes
    }

    /// Record ranges within [`bytes`](Self::bytes).
    pub fn record_refs(&self) -> &'a [FastaRecordRef] {
        self.records
    }

    /// Zero-based index of the first record in this batch.
    pub fn first_record_index(&self) -> u64 {
        self.first_record_index
    }

    /// Iterate borrowed records in this batch.
    pub fn records(&self) -> impl Iterator<Item = FastaRecord<'_>> {
        self.records
            .iter()
            .zip(self.id_tokens.iter().cloned())
            .map(|(record, id_token)| FastaRecord {
                bytes: self.bytes,
                record,
                id_token,
            })
    }

    /// Copy this borrowed batch into an owned transferable batch.
    pub fn to_owned_batch(&self) -> OwnedFastaBatch {
        OwnedFastaBatch {
            bytes: self.bytes.to_vec(),
            records: self.records.to_vec(),
            id_tokens: self.id_tokens.to_vec(),
            first_record_index: self.first_record_index,
        }
    }
}

/// Owned FASTA batch that can be moved to worker threads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OwnedFastaBatch {
    bytes: Vec<u8>,
    records: Vec<FastaRecordRef>,
    id_tokens: Vec<ByteRange>,
    first_record_index: u64,
}

impl OwnedFastaBatch {
    /// Number of records in the batch.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Whether the batch has no records.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Raw bytes backing this batch.
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// Record ranges within [`bytes`](Self::bytes).
    pub fn record_refs(&self) -> &[FastaRecordRef] {
        &self.records
    }

    /// Zero-based index of the first record in this batch.
    pub fn first_record_index(&self) -> u64 {
        self.first_record_index
    }

    /// Iterate borrowed records in this owned batch.
    pub fn records(&self) -> impl Iterator<Item = OwnedFastaRecord<'_>> {
        self.records
            .iter()
            .zip(self.id_tokens.iter().cloned())
            .map(|(record, id_token)| OwnedFastaRecord {
                bytes: &self.bytes,
                record,
                id_token,
            })
    }
}

/// Borrowed record view over an [`OwnedFastaBatch`].
#[derive(Debug, Clone, Copy)]
pub struct OwnedFastaRecord<'a> {
    bytes: &'a [u8],
    record: &'a FastaRecordRef,
    id_token: ByteRange,
}

impl<'a> OwnedFastaRecord<'a> {
    /// Return the header line including the leading `>`.
    pub fn name(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.name.clone())]
    }

    /// Return the header line without a leading `>`.
    pub fn name_without_gt(self) -> &'a [u8] {
        let name = self.name();
        name.strip_prefix(b">").unwrap_or(name)
    }

    /// Return the first whitespace-delimited identifier token.
    pub fn id_token(self) -> &'a [u8] {
        &self.bytes[self.id_token.to_usize()]
    }

    /// Return the concatenated sequence bytes.
    pub fn seq(self) -> &'a [u8] {
        &self.bytes[to_usize(self.record.seq.clone())]
    }
}

/// Slab-style streaming FASTA reader over any [`Read`] input.
///
/// The reader accepts ordinary multiline FASTA. Header lines must begin with
/// `>`, sequence lines are concatenated without newline bytes, and blank lines
/// before the first record or between records are ignored.
#[derive(Debug)]
pub struct FastaReader<R> {
    reader: BufReader<R>,
    config: FastaConfig,
    bytes: Vec<u8>,
    records: Vec<FastaRecordRef>,
    id_tokens: Vec<ByteRange>,
    line: Vec<u8>,
    pending_header: Option<PendingHeader>,
    byte_offset: u64,
    record_index: u64,
    eof: bool,
}

#[derive(Debug)]
struct PendingHeader {
    bytes: Vec<u8>,
    byte_offset: u64,
}

impl<R: Read> FastaReader<R> {
    /// Create a reader with [`FastaConfig::default`].
    pub fn new(reader: R) -> Self {
        Self::with_config(reader, FastaConfig::default())
    }

    /// Create a reader with explicit configuration.
    pub fn with_config(reader: R, config: FastaConfig) -> Self {
        let config = FastaConfig {
            batch_records: config.batch_records.max(1),
            buffer_size: config.buffer_size.max(1024),
            expected_seq_len: config.expected_seq_len,
        };
        let seq_capacity = config
            .batch_records
            .saturating_mul(config.expected_seq_len)
            .min(8 * 1024 * 1024);
        Self {
            reader: BufReader::with_capacity(config.buffer_size, reader),
            bytes: Vec::with_capacity(seq_capacity),
            records: Vec::with_capacity(config.batch_records),
            id_tokens: Vec::with_capacity(config.batch_records),
            config,
            line: Vec::new(),
            pending_header: None,
            byte_offset: 0,
            record_index: 0,
            eof: false,
        }
    }

    /// Return the wrapped reader.
    pub fn into_inner(self) -> R {
        self.reader.into_inner()
    }

    /// Read the next batch of FASTA records.
    ///
    /// Returns `Ok(None)` at EOF. A malformed non-empty line before the first
    /// header is reported as a format error.
    pub fn next_batch(&mut self) -> Result<Option<FastaBatch<'_>>> {
        self.bytes.clear();
        self.records.clear();
        self.id_tokens.clear();
        let first_record_index = self.record_index;

        while self.records.len() < self.config.batch_records {
            let Some(header) = self.next_header()? else {
                break;
            };
            self.push_record(header)?;
        }

        if self.records.is_empty() {
            return Ok(None);
        }

        self.record_index += self.records.len() as u64;
        Ok(Some(FastaBatch {
            bytes: &self.bytes,
            records: &self.records,
            id_tokens: &self.id_tokens,
            first_record_index,
        }))
    }

    /// Read the next FASTA batch into an owned transferable buffer.
    pub fn next_owned_batch(&mut self) -> Result<Option<OwnedFastaBatch>> {
        self.next_batch()
            .map(|batch| batch.map(|batch| batch.to_owned_batch()))
    }

    /// Visit every FASTA record in the stream.
    ///
    /// This is a convenience path for single-pass consumers. It uses the same
    /// batching and multiline sequence folding as [`next_batch`](Self::next_batch).
    pub fn visit_records<F>(&mut self, mut visit: F) -> Result<()>
    where
        F: FnMut(FastaVisitRecord<'_>) -> Result<()>,
    {
        self.visit_records_with_sink(&mut visit)
    }

    /// Visit every FASTA record using a sink implementation.
    pub fn visit_records_with_sink<S>(&mut self, sink: &mut S) -> Result<()>
    where
        S: FastaRecordSink,
    {
        while let Some(batch) = self.next_batch()? {
            for record in batch.records() {
                sink.record(FastaVisitRecord {
                    name: record.name(),
                    seq: record.seq(),
                })?;
            }
        }
        Ok(())
    }

    /// Count all records and bases in an ordinary FASTA stream.
    ///
    /// This uses the robust multiline FASTA parser, so it accepts wrapped
    /// sequence records and blank lines in the same way as
    /// [`next_batch`](Self::next_batch).
    pub fn stats(&mut self) -> Result<FastaStats> {
        let mut stats = FastaStats::default();
        while let Some(header) = self.next_header()? {
            self.observe_record_stats(header, &mut stats)?;
            self.record_index += 1;
        }
        Ok(stats)
    }

    fn observe_record_stats(
        &mut self,
        _header: PendingHeader,
        stats: &mut FastaStats,
    ) -> Result<()> {
        let record_index = self.record_index;
        let mut bases = 0_u64;
        let mut first = 0_u8;
        let mut last = 0_u8;

        loop {
            let line_start = self.byte_offset;
            let n = self.read_line()?;
            if n == 0 {
                self.eof = true;
                break;
            }
            let trimmed = trim_line(&self.line);
            if trimmed.starts_with(b">") {
                validate_header(trimmed, line_start, record_index + 1)?;
                self.pending_header = Some(PendingHeader {
                    bytes: trimmed.to_vec(),
                    byte_offset: line_start,
                });
                break;
            }
            if trimmed.is_empty() {
                continue;
            }
            if bases == 0 {
                first = trimmed[0];
            }
            last = *trimmed.last().unwrap_or(&0);
            bases += trimmed.len() as u64;
        }

        stats.observe_sequence_parts(bases, first, last);
        Ok(())
    }

    fn next_header(&mut self) -> Result<Option<PendingHeader>> {
        if let Some(header) = self.pending_header.take() {
            return Ok(Some(header));
        }
        if self.eof {
            return Ok(None);
        }

        loop {
            let line_start = self.byte_offset;
            let n = self.read_line()?;
            if n == 0 {
                self.eof = true;
                return Ok(None);
            }
            let trimmed = trim_line(&self.line);
            if trimmed.is_empty() {
                continue;
            }
            if !trimmed.starts_with(b">") {
                return Err(format_at(
                    "FASTA record header must start with `>`",
                    line_start,
                    self.record_index,
                ));
            }
            validate_header(trimmed, line_start, self.record_index)?;
            return Ok(Some(PendingHeader {
                bytes: trimmed.to_vec(),
                byte_offset: line_start,
            }));
        }
    }

    fn push_record(&mut self, header: PendingHeader) -> Result<()> {
        let record_index = self.record_index + self.records.len() as u64;
        let name_start = checked_u32(self.bytes.len())?;
        self.bytes.extend_from_slice(&header.bytes);
        let name_end = checked_u32(self.bytes.len())?;
        let id_token = id_token_range(&self.bytes, name_start..name_end)?;
        let seq_start = checked_u32(self.bytes.len())?;

        loop {
            let line_start = self.byte_offset;
            let n = self.read_line()?;
            if n == 0 {
                self.eof = true;
                break;
            }
            let trimmed = trim_line(&self.line);
            if trimmed.starts_with(b">") {
                validate_header(trimmed, line_start, record_index + 1)?;
                self.pending_header = Some(PendingHeader {
                    bytes: trimmed.to_vec(),
                    byte_offset: line_start,
                });
                break;
            }
            if trimmed.is_empty() {
                continue;
            }
            self.bytes.extend_from_slice(trimmed);
        }

        let seq_end = checked_u32(self.bytes.len())?;
        self.records.push(FastaRecordRef {
            name: name_start..name_end,
            seq: seq_start..seq_end,
        });
        self.id_tokens.push(id_token);
        let _ = header.byte_offset;
        Ok(())
    }

    fn read_line(&mut self) -> Result<usize> {
        self.line.clear();
        let n = self.reader.read_until(b'\n', &mut self.line)?;
        self.byte_offset += n as u64;
        Ok(n)
    }
}

/// Count records and bases from an ordinary FASTA stream.
///
/// Unlike [`count_two_line_fasta_read`], this accepts wrapped/multiline FASTA
/// using the robust [`FastaReader`] parser.
pub fn count_fasta_read<R: Read>(reader: R) -> Result<FastaStats> {
    let mut reader = FastaReader::new(reader);
    reader.stats()
}

/// Count records and bases from resident FASTA bytes.
///
/// This accepts ordinary wrapped/multiline FASTA and shares validation behavior
/// with [`visit_fasta_bytes`].
pub fn count_fasta_bytes(bytes: &[u8]) -> Result<FastaStats> {
    let mut stats = FastaStats::default();
    visit_fasta_bytes(bytes, |record| {
        stats.observe_sequence(record.seq());
        Ok(())
    })?;
    Ok(stats)
}

/// Build a `.fai`-style index over an ordinary FASTA stream.
///
/// The resulting offsets are byte offsets in the uncompressed FASTA stream.
/// Sequence records must use consistent wrapping: every non-final sequence line
/// for a record must have the same base count and byte width.
pub fn build_fasta_index<R: Read>(reader: R) -> Result<FastaIndex> {
    let mut builder = FastaIndexBuilder::default();
    let mut reader = BufReader::new(reader);
    build_fasta_index_bufread(&mut reader, &mut builder)?;
    Ok(builder.finish())
}

/// Build a `.fai`-style index over a complete BGZF-compressed FASTA stream.
///
/// The standard `.fai` offsets remain uncompressed byte offsets. Each entry also
/// carries a BGZF virtual offset for the first sequence byte, allowing callers
/// to pair the index with [`crate::BgzfSeekReader`].
#[cfg(feature = "bgzf")]
pub fn build_fasta_index_bgzf<R: Read>(reader: R) -> Result<FastaIndex> {
    let mut block_reader = BgzfDecodedBlockReader::new(reader);
    let mut builder = FastaIndexBuilder::default();
    let mut bgzf_entries = Vec::new();
    let mut line = Vec::new();

    while let Some(block) = block_reader.next_block()? {
        bgzf_entries.push(block.index_entry());
        for &byte in block.bytes() {
            line.push(byte);
            if byte == b'\n' {
                builder.observe_physical_line(&line)?;
                line.clear();
            }
        }
    }
    if !line.is_empty() {
        builder.observe_physical_line(&line)?;
    }
    builder.finish_current()?;

    let mut index = builder.finish();
    for entry in &mut index.entries {
        entry.virtual_offset = bgzf_virtual_offset_for(&bgzf_entries, entry.offset)?;
    }
    Ok(index)
}

#[cfg(feature = "bgzf")]
fn bgzf_virtual_offset_for(
    entries: &[BgzfIndexEntry],
    offset: u64,
) -> Result<Option<BgzfVirtualOffset>> {
    let idx = entries.partition_point(|entry| entry.uncompressed_offset <= offset);
    let Some(entry) = idx.checked_sub(1).and_then(|idx| entries.get(idx)) else {
        return Ok(None);
    };
    entry.virtual_offset_for(offset)
}

#[derive(Debug, Default)]
struct FastaIndexBuilder {
    entries: Vec<FastaIndexEntry>,
    seen_names: HashSet<Vec<u8>>,
    current: Option<FastaIndexRecord>,
    byte_offset: u64,
}

#[derive(Debug)]
struct FastaIndexRecord {
    name: Vec<u8>,
    len: u64,
    offset: Option<u64>,
    line_bases: Option<u64>,
    line_width: Option<u64>,
    last_line_bases: Option<u64>,
    last_line_width: Option<u64>,
    record_index: u64,
}

fn build_fasta_index_bufread<R: BufRead>(
    reader: &mut R,
    builder: &mut FastaIndexBuilder,
) -> Result<()> {
    let mut line = Vec::new();
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        builder.observe_physical_line(&line)?;
    }
    builder.finish_current()?;
    Ok(())
}

impl FastaIndexBuilder {
    fn observe_physical_line(&mut self, line: &[u8]) -> Result<()> {
        let line_start = self.byte_offset;
        self.byte_offset += line.len() as u64;
        let trimmed = trim_line(line);
        if trimmed.is_empty() {
            self.observe_blank(line_start)
        } else if trimmed.starts_with(b">") {
            self.start_record(trimmed, line_start)
        } else {
            self.observe_sequence_line(trimmed.len() as u64, line.len() as u64, line_start)
        }
    }

    fn start_record(&mut self, header: &[u8], byte_offset: u64) -> Result<()> {
        self.finish_current()?;
        validate_header(header, byte_offset, self.entries.len() as u64)?;
        let name = fasta_index_name(header);
        if self.seen_names.contains(name) {
            return Err(format_at(
                "duplicate FASTA index reference name",
                byte_offset,
                self.entries.len() as u64,
            ));
        }
        self.seen_names.insert(name.to_vec());
        self.current = Some(FastaIndexRecord {
            name: name.to_vec(),
            len: 0,
            offset: None,
            line_bases: None,
            line_width: None,
            last_line_bases: None,
            last_line_width: None,
            record_index: self.entries.len() as u64,
        });
        Ok(())
    }

    fn observe_blank(&mut self, byte_offset: u64) -> Result<()> {
        if self
            .current
            .as_ref()
            .and_then(|record| record.offset)
            .is_some()
        {
            return Err(format_at(
                "FASTA index does not support blank sequence lines",
                byte_offset,
                self.current
                    .as_ref()
                    .map_or(0, |record| record.record_index),
            ));
        }
        Ok(())
    }

    fn observe_sequence_line(&mut self, bases: u64, width: u64, byte_offset: u64) -> Result<()> {
        let Some(record) = self.current.as_mut() else {
            return Err(format_at(
                "FASTA record header must start with `>`",
                byte_offset,
                self.entries.len() as u64,
            ));
        };
        if record.offset.is_none() {
            record.offset = Some(byte_offset);
            record.line_bases = Some(bases);
            record.line_width = Some(width);
        } else if let (Some(last_bases), Some(last_width)) =
            (record.last_line_bases, record.last_line_width)
        {
            let expected_bases = record.line_bases.unwrap_or(last_bases);
            let expected_width = record.line_width.unwrap_or(last_width);
            if bases > expected_bases {
                return Err(format_at(
                    "FASTA final sequence line is longer than the first sequence line",
                    byte_offset,
                    record.record_index,
                ));
            }
            if last_bases != expected_bases || last_width != expected_width {
                return Err(format_at(
                    "non-final FASTA sequence line has inconsistent wrapping",
                    byte_offset,
                    record.record_index,
                ));
            }
        }
        record.len += bases;
        record.last_line_bases = Some(bases);
        record.last_line_width = Some(width);
        Ok(())
    }

    fn finish_current(&mut self) -> Result<()> {
        let Some(record) = self.current.take() else {
            return Ok(());
        };
        let offset = record.offset.unwrap_or(self.byte_offset);
        let line_bases = record.line_bases.unwrap_or(0);
        let line_width = record.line_width.unwrap_or(0);
        self.entries.push(FastaIndexEntry {
            name: record.name,
            len: record.len,
            offset,
            line_bases,
            line_width,
            #[cfg(feature = "bgzf")]
            virtual_offset: None,
        });
        Ok(())
    }

    fn finish(mut self) -> FastaIndex {
        let name_to_index = self
            .entries
            .iter()
            .enumerate()
            .map(|(idx, entry)| (entry.name.clone(), idx))
            .collect();
        FastaIndex {
            entries: std::mem::take(&mut self.entries),
            name_to_index,
        }
    }
}

fn fasta_index_name(header: &[u8]) -> &[u8] {
    let name = header.strip_prefix(b">").unwrap_or(header);
    let end = name
        .iter()
        .position(u8::is_ascii_whitespace)
        .unwrap_or(name.len());
    &name[..end]
}

/// Visit records from an already resident FASTA byte slice.
///
/// This path is intended for memory-mapped files, cached datasets, and other
/// callers that already own a complete FASTA byte buffer. Single-line
/// sequences are borrowed directly from the input. Multiline sequences are
/// folded into one reusable scratch buffer before the visitor is called.
///
/// Returns the number of visited records.
pub fn visit_fasta_bytes<F>(bytes: &[u8], mut visit: F) -> Result<u64>
where
    F: FnMut(FastaVisitRecord<'_>) -> Result<()>,
{
    let mut cursor = 0;
    let mut pending_header = None;
    let mut record_index = 0;
    let mut folded = Vec::new();

    while let Some((header_offset, header)) =
        take_next_header(bytes, &mut cursor, &mut pending_header, record_index)?
    {
        validate_header(header, header_offset, record_index)?;

        folded.clear();
        let mut first_seq = None;

        loop {
            if cursor >= bytes.len() {
                break;
            }
            let line_start = cursor as u64;
            let line = next_trimmed_line(bytes, &mut cursor);
            if line.starts_with(b">") {
                validate_header(line, line_start, record_index + 1)?;
                pending_header = Some((line_start, line));
                break;
            }
            if line.is_empty() {
                continue;
            }

            first_seq = Some(line);
            break;
        }

        let Some(first_seq) = first_seq else {
            visit(FastaVisitRecord {
                name: header,
                seq: b"",
            })?;
            record_index += 1;
            continue;
        };

        if cursor >= bytes.len() || bytes[cursor] == b'>' {
            visit(FastaVisitRecord {
                name: header,
                seq: first_seq,
            })?;
            record_index += 1;
            continue;
        }

        let mut seq_line_count = 1_usize;
        loop {
            if cursor >= bytes.len() {
                break;
            }
            let line_start = cursor as u64;
            let line = next_trimmed_line(bytes, &mut cursor);
            if line.starts_with(b">") {
                validate_header(line, line_start, record_index + 1)?;
                pending_header = Some((line_start, line));
                break;
            }
            if line.is_empty() {
                continue;
            }

            seq_line_count += 1;
            if seq_line_count == 2 {
                folded.extend_from_slice(first_seq);
            }
            folded.extend_from_slice(line);
        }

        let seq = if seq_line_count == 1 {
            first_seq
        } else {
            &folded
        };
        visit(FastaVisitRecord { name: header, seq })?;
        record_index += 1;
    }

    Ok(record_index)
}

/// Detect whether resident FASTA bytes are strict two-line FASTA.
///
/// The detector validates headers enough to reject non-FASTA leading content
/// and empty headers. It returns [`FastaShape::Multiline`] for valid FASTA that
/// needs the robust parser, including blank lines between records.
pub fn detect_fasta_shape(bytes: &[u8]) -> Result<FastaShape> {
    let mut cursor = 0;
    let mut saw_record = false;

    while cursor < bytes.len() {
        let header_offset = cursor as u64;
        let header = next_trimmed_line(bytes, &mut cursor);
        if header.is_empty() {
            if saw_record {
                return Ok(FastaShape::Multiline);
            }
            continue;
        }
        if !header.starts_with(b">") {
            return Err(format_at(
                "FASTA record header must start with `>`",
                header_offset,
                0,
            ));
        }
        validate_header(header, header_offset, 0)?;
        saw_record = true;

        let mut seq_lines = 0_usize;
        while cursor < bytes.len() {
            let checkpoint = cursor;
            let line = next_trimmed_line(bytes, &mut cursor);
            if line.starts_with(b">") {
                cursor = checkpoint;
                break;
            }
            if line.is_empty() {
                return Ok(FastaShape::Multiline);
            }
            seq_lines += 1;
            if seq_lines > 1 {
                return Ok(FastaShape::Multiline);
            }
        }
        if seq_lines != 1 {
            return Ok(FastaShape::Multiline);
        }
    }

    if saw_record {
        Ok(FastaShape::TwoLine)
    } else {
        Ok(FastaShape::Empty)
    }
}

/// Visit resident FASTA bytes with automatic two-line fast-path detection.
pub fn visit_fasta_bytes_auto<F>(bytes: &[u8], visit: F) -> Result<u64>
where
    F: FnMut(FastaVisitRecord<'_>) -> Result<()>,
{
    match detect_fasta_shape(bytes)? {
        FastaShape::TwoLine => visit_two_line_fasta_bytes(bytes, visit),
        FastaShape::Empty | FastaShape::Multiline => visit_fasta_bytes(bytes, visit),
    }
}

/// Count strict two-line FASTA records from resident bytes.
///
/// This is the lightest resident path for canonical `>header\nsequence\n`
/// FASTA when callers only need record counts, total bases, and a deterministic
/// shape checksum. It validates the same strict two-line structure as
/// [`visit_two_line_fasta_bytes`].
pub fn count_two_line_fasta_bytes(bytes: &[u8]) -> Result<FastaStats> {
    let mut cursor = 0;
    let mut record_index = 0;
    let mut stats = FastaStats::default();

    while cursor < bytes.len() {
        let header_offset = cursor as u64;
        if bytes[cursor] != b'>' {
            return Err(format_at(
                "two-line FASTA record header must start with `>`",
                header_offset,
                record_index,
            ));
        }
        let header = next_trimmed_line(bytes, &mut cursor);
        validate_header(header, header_offset, record_index)?;
        if cursor >= bytes.len() {
            return Err(format_at(
                "two-line FASTA record is missing a sequence line",
                cursor as u64,
                record_index,
            ));
        }

        let seq_offset = cursor as u64;
        let seq = next_trimmed_line(bytes, &mut cursor);
        if seq.is_empty() || seq.starts_with(b">") {
            return Err(format_at(
                "two-line FASTA record is missing a sequence line",
                seq_offset,
                record_index,
            ));
        }
        if cursor < bytes.len() && bytes[cursor] != b'>' {
            return Err(format_at(
                "two-line FASTA sequence must be followed by a header",
                cursor as u64,
                record_index + 1,
            ));
        }

        stats.observe_sequence(seq);
        record_index += 1;
    }

    Ok(stats)
}

/// Visit records from a resident, strict two-line FASTA byte slice.
///
/// This is the fastest resident FASTA path for canonical files shaped as
/// `>header\nsequence\n` repeated. It rejects blank lines and multiline
/// sequence records. Use [`visit_fasta_bytes`] when ordinary multiline FASTA
/// support is required.
///
/// Returns the number of visited records.
pub fn visit_two_line_fasta_bytes<F>(bytes: &[u8], mut visit: F) -> Result<u64>
where
    F: FnMut(FastaVisitRecord<'_>) -> Result<()>,
{
    let mut cursor = 0;
    let mut record_index = 0;

    while cursor < bytes.len() {
        let header_offset = cursor as u64;
        if bytes[cursor] != b'>' {
            return Err(format_at(
                "two-line FASTA record header must start with `>`",
                header_offset,
                record_index,
            ));
        }
        let header = next_trimmed_line(bytes, &mut cursor);
        validate_header(header, header_offset, record_index)?;
        if cursor >= bytes.len() {
            return Err(format_at(
                "two-line FASTA record is missing a sequence line",
                cursor as u64,
                record_index,
            ));
        }

        let seq_offset = cursor as u64;
        let seq = next_trimmed_line(bytes, &mut cursor);
        if seq.is_empty() || seq.starts_with(b">") {
            return Err(format_at(
                "two-line FASTA record is missing a sequence line",
                seq_offset,
                record_index,
            ));
        }
        if cursor < bytes.len() && bytes[cursor] != b'>' {
            return Err(format_at(
                "two-line FASTA sequence must be followed by a header",
                cursor as u64,
                record_index + 1,
            ));
        }

        visit(FastaVisitRecord { name: header, seq })?;
        record_index += 1;
    }

    Ok(record_index)
}

/// Visit records from a strict two-line FASTA stream.
///
/// This path parses canonical `>header\nsequence\n` FASTA directly from a
/// buffered stream. Complete sequence lines are borrowed from the read buffer;
/// only headers and chunk-boundary fragments are copied. It rejects blank lines
/// and multiline sequence records. Use [`FastaReader`] or [`visit_fasta_bytes`]
/// when ordinary multiline FASTA support is required.
///
/// Returns the number of visited records.
pub fn visit_two_line_fasta_read<R, F>(reader: R, visit: F) -> Result<u64>
where
    R: Read,
    F: FnMut(FastaVisitRecord<'_>) -> Result<()>,
{
    let mut reader = BufReader::with_capacity(TWO_LINE_STREAM_BUFFER_SIZE, reader);
    visit_two_line_fasta_bufread(&mut reader, visit)
}

/// Count strict two-line FASTA records from a stream.
///
/// This is the lightest path for workloads that only need record counts, total
/// bases, and a deterministic shape checksum. It validates the canonical
/// `>header\nsequence\n` structure and rejects multiline FASTA.
pub fn count_two_line_fasta_read<R: Read>(reader: R) -> Result<FastaStats> {
    let mut reader = BufReader::with_capacity(TWO_LINE_STREAM_BUFFER_SIZE, reader);
    count_two_line_fasta_bufread(&mut reader)
}

fn visit_two_line_fasta_bufread<R, F>(reader: &mut R, mut visit: F) -> Result<u64>
where
    R: BufRead,
    F: FnMut(FastaVisitRecord<'_>) -> Result<()>,
{
    let mut state = TwoLineStreamState::Header;
    let mut record_index = 0;
    let mut header = Vec::new();
    let mut carry = Vec::new();

    loop {
        let available = reader.fill_buf()?;
        if available.is_empty() {
            break;
        }

        let mut consumed = 0;
        while consumed < available.len() {
            let Some(relative_newline) = memchr(b'\n', &available[consumed..]) else {
                carry.extend_from_slice(&available[consumed..]);
                consumed = available.len();
                break;
            };
            let line_end = consumed + relative_newline;
            let line = &available[consumed..line_end];
            process_two_line_stream_line(
                line,
                &mut carry,
                &mut state,
                &mut header,
                &mut record_index,
                &mut visit,
            )?;
            consumed = line_end + 1;
        }
        reader.consume(consumed);
    }

    if !carry.is_empty() {
        process_two_line_stream_line(
            b"",
            &mut carry,
            &mut state,
            &mut header,
            &mut record_index,
            &mut visit,
        )?;
    }

    match state {
        TwoLineStreamState::Header => Ok(record_index),
        TwoLineStreamState::Seq => Err(format_at(
            "two-line FASTA record is missing a sequence line",
            0,
            record_index,
        )),
    }
}

fn count_two_line_fasta_bufread<R: BufRead>(reader: &mut R) -> Result<FastaStats> {
    let mut stats = FastaStats::default();
    let mut state = TwoLineStreamState::Header;
    let mut record_index = 0;
    let mut carry = Vec::new();

    loop {
        let available = reader.fill_buf()?;
        if available.is_empty() {
            break;
        }

        let mut consumed = 0;
        while consumed < available.len() {
            let Some(relative_newline) = memchr(b'\n', &available[consumed..]) else {
                carry.extend_from_slice(&available[consumed..]);
                consumed = available.len();
                break;
            };
            let line_end = consumed + relative_newline;
            let line = &available[consumed..line_end];
            process_two_line_count_line(
                line,
                &mut carry,
                &mut state,
                &mut record_index,
                &mut stats,
            )?;
            consumed = line_end + 1;
        }
        reader.consume(consumed);
    }

    if !carry.is_empty() {
        process_two_line_count_line(b"", &mut carry, &mut state, &mut record_index, &mut stats)?;
    }

    if state == TwoLineStreamState::Seq {
        return Err(format_at(
            "two-line FASTA record is missing a sequence line",
            0,
            record_index,
        ));
    }
    Ok(stats)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TwoLineStreamState {
    Header,
    Seq,
}

fn process_two_line_stream_line<F>(
    line: &[u8],
    carry: &mut Vec<u8>,
    state: &mut TwoLineStreamState,
    header: &mut Vec<u8>,
    record_index: &mut u64,
    visit: &mut F,
) -> Result<()>
where
    F: FnMut(FastaVisitRecord<'_>) -> Result<()>,
{
    if carry.is_empty() {
        process_complete_two_line_stream_line(line, state, header, record_index, visit)
    } else {
        carry.extend_from_slice(line);
        let owned_line = trim_line(carry);
        process_complete_two_line_stream_line(owned_line, state, header, record_index, visit)?;
        carry.clear();
        Ok(())
    }
}

fn process_complete_two_line_stream_line<F>(
    line: &[u8],
    state: &mut TwoLineStreamState,
    header: &mut Vec<u8>,
    record_index: &mut u64,
    visit: &mut F,
) -> Result<()>
where
    F: FnMut(FastaVisitRecord<'_>) -> Result<()>,
{
    let line = trim_line(line);
    if *state == TwoLineStreamState::Header {
        if !line.starts_with(b">") {
            return Err(format_at(
                "two-line FASTA record header must start with `>`",
                0,
                *record_index,
            ));
        }
        validate_header(line, 0, *record_index)?;
        header.clear();
        header.extend_from_slice(line);
        *state = TwoLineStreamState::Seq;
        return Ok(());
    }

    if line.is_empty() || line.starts_with(b">") {
        return Err(format_at(
            "two-line FASTA record is missing a sequence line",
            0,
            *record_index,
        ));
    }
    visit(FastaVisitRecord {
        name: header,
        seq: line,
    })?;
    *record_index += 1;
    *state = TwoLineStreamState::Header;
    Ok(())
}

fn process_two_line_count_line(
    line: &[u8],
    carry: &mut Vec<u8>,
    state: &mut TwoLineStreamState,
    record_index: &mut u64,
    stats: &mut FastaStats,
) -> Result<()> {
    if carry.is_empty() {
        process_complete_two_line_count_line(line, state, record_index, stats)
    } else {
        carry.extend_from_slice(line);
        let owned_line = trim_line(carry);
        process_complete_two_line_count_line(owned_line, state, record_index, stats)?;
        carry.clear();
        Ok(())
    }
}

fn process_complete_two_line_count_line(
    line: &[u8],
    state: &mut TwoLineStreamState,
    record_index: &mut u64,
    stats: &mut FastaStats,
) -> Result<()> {
    let line = trim_line(line);
    if *state == TwoLineStreamState::Header {
        if !line.starts_with(b">") {
            return Err(format_at(
                "two-line FASTA record header must start with `>`",
                0,
                *record_index,
            ));
        }
        validate_header(line, 0, *record_index)?;
        *state = TwoLineStreamState::Seq;
        return Ok(());
    }

    if line.is_empty() || line.starts_with(b">") {
        return Err(format_at(
            "two-line FASTA record is missing a sequence line",
            0,
            *record_index,
        ));
    }
    stats.observe_sequence(line);
    *record_index += 1;
    *state = TwoLineStreamState::Header;
    Ok(())
}

fn take_next_header<'a>(
    bytes: &'a [u8],
    cursor: &mut usize,
    pending_header: &mut Option<(u64, &'a [u8])>,
    record_index: u64,
) -> Result<Option<(u64, &'a [u8])>> {
    if let Some(header) = pending_header.take() {
        return Ok(Some(header));
    }

    while *cursor < bytes.len() {
        let line_start = *cursor as u64;
        let line = next_trimmed_line(bytes, cursor);
        if line.is_empty() {
            continue;
        }
        if !line.starts_with(b">") {
            return Err(format_at(
                "FASTA record header must start with `>`",
                line_start,
                record_index,
            ));
        }
        return Ok(Some((line_start, line)));
    }

    Ok(None)
}

fn next_trimmed_line<'a>(bytes: &'a [u8], cursor: &mut usize) -> &'a [u8] {
    let start = *cursor;
    match memchr(b'\n', &bytes[start..]) {
        Some(relative) => {
            let end = start + relative;
            *cursor = end + 1;
            trim_line(&bytes[start..end])
        }
        None => {
            *cursor = bytes.len();
            trim_line(&bytes[start..])
        }
    }
}

fn trim_line(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    line.strip_suffix(b"\r").unwrap_or(line)
}

fn validate_header(header: &[u8], byte_offset: u64, record_index: u64) -> Result<()> {
    if header.len() == 1 || fasta_index_name(header).is_empty() {
        return Err(format_at("empty FASTA id", byte_offset, record_index));
    }
    Ok(())
}

fn id_token_range(bytes: &[u8], name: Range<u32>) -> Result<ByteRange> {
    let name_range = to_usize(name.clone());
    let header = bytes
        .get(name_range)
        .ok_or_else(|| FastqError::Format("FASTA header byte range exceeds batch buffer".into()))?;
    let without_gt = header.strip_prefix(b">").unwrap_or(header);
    let token_end = without_gt
        .iter()
        .position(u8::is_ascii_whitespace)
        .unwrap_or(without_gt.len());
    let start = name.start + u32::from(header.starts_with(b">"));
    let token_end_u32 = u32::try_from(token_end)
        .map_err(|_| FastqError::Format("FASTA id token range exceeds u32 range".into()))?;
    let end = start
        .checked_add(token_end_u32)
        .ok_or_else(|| FastqError::Format("FASTA id token range overflowed".into()))?;
    Ok(ByteRange { start, end })
}

fn checked_u32(value: usize) -> Result<u32> {
    u32::try_from(value)
        .map_err(|_| FastqError::Format("FASTA batch byte offsets exceed u32 range".into()))
}

fn format_at(message: impl Into<String>, byte_offset: u64, record_index: u64) -> FastqError {
    FastqError::FormatAt {
        message: message.into(),
        position: FastqPosition::new(byte_offset, record_index, 0),
    }
}

fn to_usize(range: Range<u32>) -> Range<usize> {
    range.start as usize..range.end as usize
}

#[cfg(test)]
mod tests;
