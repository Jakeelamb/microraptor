use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::ops::Range;

use crate::error::{FastqError, FastqPosition, Result};
#[cfg(feature = "bgzf")]
use crate::{BgzfDecodedBlockReader, BgzfIndex, BgzfIndexEntry, BgzfSeekReader, BgzfVirtualOffset};
use memchr::memchr;

const DEFAULT_BATCH_RECORDS: usize = 1024;
const DEFAULT_READER_BUFFER_SIZE: usize = 64 * 1024;
const TWO_LINE_STREAM_BUFFER_SIZE: usize = 64 * 1024;

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
        let name = self.name_without_gt();
        let end = name
            .iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(name.len());
        &name[..end]
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
        Ok(self.offset + (pos / self.line_bases) * self.line_width + (pos % self.line_bases))
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
            spans.push(start..start + take);
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

#[cfg(feature = "bgzf")]
fn copy_exact_into<R: Read>(
    reader: &mut R,
    mut len: u64,
    out: &mut Vec<u8>,
    scratch: &mut [u8],
) -> Result<()> {
    while len > 0 {
        let take = usize::try_from(len.min(scratch.len() as u64))
            .map_err(|_| FastqError::Format("FASTA fetch span exceeds usize range".into()))?;
        reader.read_exact(&mut scratch[..take])?;
        out.extend_from_slice(&scratch[..take]);
        len -= take as u64;
    }
    Ok(())
}

/// Seekable FASTA reader backed by a `.fai` index.
pub struct IndexedFastaReader<R> {
    inner: R,
    index: FastaIndex,
}

impl<R: Read + Seek> IndexedFastaReader<R> {
    /// Create a seekable FASTA reader from an input stream and index.
    pub fn new(inner: R, index: FastaIndex) -> Self {
        Self { inner, index }
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
        entry.validate_range(range.clone())?;
        out.clear();
        out.reserve(usize::try_from(range.end - range.start).map_err(|_| {
            FastqError::Format("FASTA fetch range length exceeds usize range".into())
        })?);
        for span in entry.sequence_spans(range)? {
            self.inner.seek(SeekFrom::Start(span.start))?;
            let len = usize::try_from(span.end - span.start).map_err(|_| {
                FastqError::Format("FASTA physical span length exceeds usize range".into())
            })?;
            let start = out.len();
            out.resize(start + len, 0);
            self.inner.read_exact(&mut out[start..])?;
        }
        Ok(())
    }

    /// Fetch a zero-based half-open sequence range into an owned buffer.
    pub fn fetch(&mut self, name: &[u8], range: Range<u64>) -> Result<Vec<u8>> {
        let mut out = Vec::new();
        self.fetch_into(name, range, &mut out)?;
        Ok(out)
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
}

#[cfg(feature = "bgzf")]
impl<R: Read + Seek> BgzfIndexedFastaReader<R> {
    /// Create a BGZF FASTA random-access reader.
    pub fn new(inner: R, fasta_index: FastaIndex, bgzf_index: BgzfIndex) -> Self {
        Self {
            inner: BgzfSeekReader::new(inner),
            fasta_index,
            bgzf_index,
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
        entry.validate_range(range.clone())?;
        out.clear();
        out.reserve(usize::try_from(range.end - range.start).map_err(|_| {
            FastqError::Format("FASTA fetch range length exceeds usize range".into())
        })?);
        let mut scratch = vec![0_u8; 8192];
        for span in entry.sequence_spans(range)? {
            let virtual_offset = self
                .bgzf_index
                .virtual_offset_for_uncompressed_offset(span.start)?
                .ok_or_else(|| FastqError::Bgzf("BGZF span offset is not indexed".into()))?;
            self.inner.seek_virtual_offset(virtual_offset)?;
            copy_exact_into(&mut self.inner, span.end - span.start, out, &mut scratch)?;
        }
        Ok(())
    }

    /// Fetch a zero-based half-open sequence range into an owned buffer.
    pub fn fetch(&mut self, name: &[u8], range: Range<u64>) -> Result<Vec<u8>> {
        let mut out = Vec::new();
        self.fetch_into(name, range, &mut out)?;
        Ok(out)
    }
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
        self.records.iter().map(|record| FastaRecord {
            bytes: self.bytes,
            record,
        })
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
            first_record_index,
        }))
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
    visit_two_line_fasta_bufread(reader, |record| {
        stats.observe_sequence(record.seq());
        Ok(())
    })?;
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
    if header.len() == 1 {
        return Err(format_at("empty FASTA id", byte_offset, record_index));
    }
    Ok(())
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
mod tests {
    use super::*;

    #[test]
    fn reads_multiline_fasta_records() {
        let input = b">seq1 description\nACG\nTN\n>seq2\nGG\n";
        let mut reader = FastaReader::new(&input[..]);
        let batch = reader.next_batch().unwrap().unwrap();
        let records: Vec<_> = batch.records().collect();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name(), b">seq1 description");
        assert_eq!(records[0].id_token(), b"seq1");
        assert_eq!(records[0].seq(), b"ACGTN");
        assert_eq!(records[1].name_without_gt(), b"seq2");
        assert_eq!(records[1].seq(), b"GG");
        assert!(reader.next_batch().unwrap().is_none());
    }

    #[test]
    fn carries_header_across_batches() {
        let input = b">seq1\nAC\n>seq2\nGT\n";
        let mut reader = FastaReader::with_config(
            &input[..],
            FastaConfig {
                batch_records: 1,
                ..FastaConfig::default()
            },
        );
        let first = reader.next_batch().unwrap().unwrap();
        let first_record = first.records().next().unwrap();
        assert_eq!(first_record.name_without_gt(), b"seq1");
        assert_eq!(first_record.seq(), b"AC");

        let second = reader.next_batch().unwrap().unwrap();
        let second_record = second.records().next().unwrap();
        assert_eq!(second_record.name_without_gt(), b"seq2");
        assert_eq!(second_record.seq(), b"GT");
        assert_eq!(second.first_record_index(), 1);
    }

    #[test]
    fn rejects_non_header_before_first_record() {
        let mut reader = FastaReader::new(&b"ACGT\n"[..]);
        let err = reader.next_batch().unwrap_err();
        assert!(err.to_string().contains("FASTA record header"));
    }

    #[test]
    fn rejects_empty_header() {
        let mut reader = FastaReader::new(&b">\nACGT\n"[..]);
        let err = reader.next_batch().unwrap_err();
        assert!(err.to_string().contains("empty FASTA id"));
    }

    #[test]
    fn visits_records() {
        let input = b">seq1\nAC\n>seq2\nGT\n";
        let mut reader = FastaReader::new(&input[..]);
        let mut seen = Vec::new();
        reader
            .visit_records(|record| {
                seen.push((record.id_token().to_vec(), record.seq().to_vec()));
                Ok(())
            })
            .unwrap();

        assert_eq!(
            seen,
            vec![
                (b"seq1".to_vec(), b"AC".to_vec()),
                (b"seq2".to_vec(), b"GT".to_vec())
            ]
        );
    }

    #[test]
    fn reports_later_empty_header_index() {
        let input = b">seq1\nAC\n>seq2\nGT\n>\nTT\n";
        let mut reader = FastaReader::with_config(
            &input[..],
            FastaConfig {
                batch_records: 8,
                ..FastaConfig::default()
            },
        );
        let err = reader.next_batch().unwrap_err();
        assert!(err.to_string().contains("record 2"));
    }

    #[test]
    fn visits_resident_fasta_bytes() {
        let input = b"\n>seq1 description\nAC\nGT\n>seq2\nTT\n";
        let mut seen = Vec::new();
        let records = visit_fasta_bytes(input, |record| {
            seen.push((record.id_token().to_vec(), record.seq().to_vec()));
            Ok(())
        })
        .unwrap();

        assert_eq!(records, 2);
        assert_eq!(
            seen,
            vec![
                (b"seq1".to_vec(), b"ACGT".to_vec()),
                (b"seq2".to_vec(), b"TT".to_vec())
            ]
        );
    }

    #[test]
    fn resident_visitor_rejects_non_header_before_first_record() {
        let err = visit_fasta_bytes(b"ACGT\n", |_| Ok(())).unwrap_err();
        assert!(err.to_string().contains("FASTA record header"));
    }

    #[test]
    fn resident_visitor_reports_later_empty_header_index() {
        let input = b">seq1\nAC\n>\nTT\n";
        let err = visit_fasta_bytes(input, |_| Ok(())).unwrap_err();
        assert!(err.to_string().contains("record 1"));
    }

    #[test]
    fn detects_fasta_shape() {
        assert_eq!(detect_fasta_shape(b"").unwrap(), FastaShape::Empty);
        assert_eq!(
            detect_fasta_shape(b">seq1\nAC\n>seq2\nTT\n").unwrap(),
            FastaShape::TwoLine
        );
        assert_eq!(
            detect_fasta_shape(b">seq1\nAC\nGT\n").unwrap(),
            FastaShape::Multiline
        );
        assert_eq!(
            detect_fasta_shape(b">seq1\nAC\n\n>seq2\nTT\n").unwrap(),
            FastaShape::Multiline
        );
        assert!(detect_fasta_shape(b"ACGT\n").is_err());
    }

    #[test]
    fn auto_resident_visitor_matches_robust_multiline_visitor() {
        let input = b">seq1\nAC\nGT\n>seq2\nTT\n";
        let mut auto = Vec::new();
        let mut robust = Vec::new();

        visit_fasta_bytes_auto(input, |record| {
            auto.push((record.id_token().to_vec(), record.seq().to_vec()));
            Ok(())
        })
        .unwrap();
        visit_fasta_bytes(input, |record| {
            robust.push((record.id_token().to_vec(), record.seq().to_vec()));
            Ok(())
        })
        .unwrap();

        assert_eq!(auto, robust);
    }

    #[test]
    fn counts_wrapped_fasta_read_and_bytes() {
        let input = b">seq1\nAC\nGT\n>seq2\nTTA\nA\n";
        let from_read = count_fasta_read(&input[..]).unwrap();
        let from_bytes = count_fasta_bytes(input).unwrap();
        let mut reader = FastaReader::new(&input[..]);
        let from_reader = reader.stats().unwrap();

        assert_eq!(from_read, from_bytes);
        assert_eq!(from_reader, from_read);
        assert_eq!(from_read.records, 2);
        assert_eq!(from_read.bases, 8);
    }

    #[test]
    fn stats_continue_after_partial_batch_without_folding_next_records() {
        let input = b">seq1\nAC\n>seq2\nGT\nTA\n>seq3\nCC\n";
        let mut reader = FastaReader::with_config(
            &input[..],
            FastaConfig {
                batch_records: 1,
                ..FastaConfig::default()
            },
        );
        {
            let first = reader.next_batch().unwrap().unwrap();
            assert_eq!(first.records().next().unwrap().seq(), b"AC");
        }

        let stats = reader.stats().unwrap();
        let mut expected = FastaStats::default();
        expected.observe_sequence(b"GTTA");
        expected.observe_sequence(b"CC");
        assert_eq!(stats, expected);
    }

    #[test]
    fn builds_fasta_index_for_wrapped_reference() {
        let input = b">chr1 description\nACGT\nAC\n>chr2\nTTTT\n";
        let index = build_fasta_index(&input[..]).unwrap();
        assert_eq!(index.len(), 2);

        let chr1 = index.get(b"chr1").unwrap();
        assert_eq!(chr1.len, 6);
        assert_eq!(chr1.offset, 18);
        assert_eq!(chr1.line_bases, 4);
        assert_eq!(chr1.line_width, 5);

        let chr2 = index.get(b"chr2").unwrap();
        assert_eq!(chr2.len, 4);
        assert_eq!(chr2.line_bases, 4);
        assert_eq!(chr2.line_width, 5);

        assert_eq!(
            index.to_fai_string(),
            "chr1\t6\t18\t4\t5\nchr2\t4\t32\t4\t5\n"
        );
    }

    #[test]
    fn fasta_index_rejects_inconsistent_non_final_wrapping() {
        let err = build_fasta_index(&b">chr1\nAC\nACGT\nA\n"[..]).unwrap_err();
        assert!(err.to_string().contains("longer than the first"));
    }

    #[test]
    fn fasta_index_rejects_short_internal_wrapping() {
        let err = build_fasta_index(&b">chr1\nACGT\nAC\nA\n"[..]).unwrap_err();
        assert!(err.to_string().contains("inconsistent wrapping"));
    }

    #[test]
    fn fasta_index_rejects_duplicate_names() {
        let err = build_fasta_index(&b">chr1\nAC\n>chr1 desc\nGT\n"[..]).unwrap_err();
        assert!(err.to_string().contains("duplicate"));
    }

    #[test]
    #[cfg(feature = "bgzf")]
    fn builds_bgzf_aware_fasta_index() {
        let input = b">chr1\nACGT\nAC\n>chr2\nTTTT\n";
        let encoded = crate::compress_bgzf_parallel(input, 2).unwrap();
        let index = build_fasta_index_bgzf(&encoded[..]).unwrap();
        let chr1 = index.get(b"chr1").unwrap();
        assert_eq!(chr1.len, 6);
        let vo = chr1.virtual_offset.unwrap();
        assert_eq!(vo.compressed_offset(), 0);
        assert_eq!(vo.in_block_offset(), 6);
    }

    #[test]
    #[cfg(feature = "bgzf")]
    fn bgzf_fasta_index_streams_lines_across_block_boundaries() {
        use std::io::Read;

        let mut input = b">chr1\n".to_vec();
        input.extend(std::iter::repeat_n(b'A', 70_000));
        input.extend_from_slice(b"\n>chr2\nTTTT\n");
        let encoded = crate::compress_bgzf_parallel(&input, 2).unwrap();
        let index = build_fasta_index_bgzf(&encoded[..]).unwrap();

        let chr1 = index.get(b"chr1").unwrap();
        assert_eq!(chr1.len, 70_000);
        assert_eq!(chr1.offset, 6);
        assert_eq!(chr1.line_bases, 70_000);
        assert_eq!(chr1.line_width, 70_001);

        let chr2 = index.get(b"chr2").unwrap();
        assert_eq!(chr2.len, 4);
        let chr2_vo = chr2.virtual_offset.unwrap();
        assert!(chr2_vo.compressed_offset() > 0);

        let mut reader = crate::BgzfSeekReader::new(std::io::Cursor::new(encoded));
        reader.seek_virtual_offset(chr2_vo).unwrap();
        let mut out = [0_u8; 4];
        reader.read_exact(&mut out).unwrap();
        assert_eq!(&out, b"TTTT");
    }

    #[test]
    fn visits_two_line_resident_fasta_bytes() {
        let input = b">seq1 description\nACGT\n>seq2\nTT\n";
        let mut seen = Vec::new();
        let records = visit_two_line_fasta_bytes(input, |record| {
            seen.push((record.id_token().to_vec(), record.seq().to_vec()));
            Ok(())
        })
        .unwrap();

        assert_eq!(records, 2);
        assert_eq!(
            seen,
            vec![
                (b"seq1".to_vec(), b"ACGT".to_vec()),
                (b"seq2".to_vec(), b"TT".to_vec())
            ]
        );
    }

    #[test]
    fn counts_two_line_resident_fasta_bytes() {
        let input = b">seq1\nACGT\n>seq2\nTT\n";
        let stats = count_two_line_fasta_bytes(input).unwrap();
        let mut reference = FastaStats::default();
        visit_two_line_fasta_bytes(input, |record| {
            reference.observe_sequence(record.seq());
            Ok(())
        })
        .unwrap();

        assert_eq!(stats, reference);
        assert_eq!(stats.records, 2);
        assert_eq!(stats.bases, 6);
    }

    #[test]
    fn two_line_visitor_rejects_multiline_fasta() {
        let err = visit_two_line_fasta_bytes(b">seq1\nAC\nGT\n", |_| Ok(())).unwrap_err();
        assert!(err.to_string().contains("followed by a header"));
    }

    #[test]
    fn two_line_resident_counter_rejects_multiline_fasta() {
        let err = count_two_line_fasta_bytes(b">seq1\nAC\nGT\n").unwrap_err();
        assert!(err.to_string().contains("followed by a header"));
    }

    #[test]
    fn visits_two_line_fasta_stream() {
        let input = b">seq1 description\nACGT\n>seq2\nTT\n";
        let mut seen = Vec::new();
        let records = visit_two_line_fasta_read(&input[..], |record| {
            seen.push((record.id_token().to_vec(), record.seq().to_vec()));
            Ok(())
        })
        .unwrap();

        assert_eq!(records, 2);
        assert_eq!(
            seen,
            vec![
                (b"seq1".to_vec(), b"ACGT".to_vec()),
                (b"seq2".to_vec(), b"TT".to_vec())
            ]
        );
    }

    #[test]
    fn counts_two_line_fasta_stream() {
        let input = b">seq1\nACGT\n>seq2\nTT\n";
        let stats = count_two_line_fasta_read(&input[..]).unwrap();
        let mut reference = FastaStats::default();
        visit_two_line_fasta_read(&input[..], |record| {
            reference.observe_sequence(record.seq());
            Ok(())
        })
        .unwrap();

        assert_eq!(stats, reference);
        assert_eq!(stats.records, 2);
        assert_eq!(stats.bases, 6);
    }

    #[test]
    fn two_line_counter_rejects_multiline_fasta() {
        let err = count_two_line_fasta_read(&b">seq1\nAC\nGT\n"[..]).unwrap_err();
        assert!(err.to_string().contains("header must start"));
    }

    #[test]
    fn two_line_stream_rejects_multiline_fasta() {
        let err = visit_two_line_fasta_read(&b">seq1\nAC\nGT\n"[..], |_| Ok(())).unwrap_err();
        assert!(err.to_string().contains("header must start"));
    }

    #[test]
    fn parses_fai_and_fetches_wrapped_range() {
        let input = b">chr1 desc\nACGT\nTGCA\nAA\n>chr2\nGG\n";
        let index = build_fasta_index(&input[..]).unwrap();
        let fai = index.to_fai_string();
        let parsed = FastaIndex::from_fai_str(&fai).unwrap();
        assert_eq!(parsed.to_fai_string(), fai);

        let chr1 = parsed.get(b"chr1").unwrap();
        assert_eq!(chr1.sequence_offset(0).unwrap(), 11);
        assert_eq!(chr1.sequence_offset(4).unwrap(), 16);
        assert_eq!(chr1.sequence_spans(2..8).unwrap(), vec![13..15, 16..20]);

        let mut reader = IndexedFastaReader::new(std::io::Cursor::new(input), parsed);
        assert_eq!(reader.fetch(b"chr1", 2..8).unwrap(), b"GTTGCA");
        assert_eq!(reader.fetch(b"chr1", 10..10).unwrap(), b"");
        assert_eq!(reader.fetch(b"chr2", 0..2).unwrap(), b"GG");
    }

    #[test]
    fn rejects_bad_fai_and_bad_fetch_ranges() {
        assert!(FastaIndex::from_fai_str("chr1\t1\t2\t3\n").is_err());
        assert!(FastaIndex::from_fai_str("chr1\t1\t2\t0\t1\n").is_err());
        assert!(FastaIndex::from_fai_str("chr1\t1\t2\t1\t1\nchr1\t1\t2\t1\t1\n").is_err());

        let input = b">chr1\nACGT\n";
        let index = build_fasta_index(&input[..]).unwrap();
        let mut reader = IndexedFastaReader::new(std::io::Cursor::new(input), index);
        assert!(reader.fetch(b"missing", 0..1).is_err());
        assert!(reader.fetch(b"chr1", Range { start: 3, end: 2 }).is_err());
        assert!(reader.fetch(b"chr1", 0..5).is_err());
    }

    #[cfg(feature = "bgzf")]
    #[test]
    fn fetches_bgzf_fasta_range_using_arbitrary_virtual_offsets() {
        let mut input = b">chr1\n".to_vec();
        let seq = (0..70_010).map(|i| b"ACGT"[i % 4]).collect::<Vec<_>>();
        for chunk in seq.chunks(80) {
            input.extend_from_slice(chunk);
            input.push(b'\n');
        }
        let encoded = crate::compress_bgzf_parallel(&input, 2).unwrap();
        let fasta_index = build_fasta_index_bgzf(&encoded[..]).unwrap();
        let bgzf_index = crate::build_bgzf_index_strict(&encoded[..]).unwrap();
        let mut reader =
            BgzfIndexedFastaReader::new(std::io::Cursor::new(encoded), fasta_index, bgzf_index);

        let fetched = reader.fetch(b"chr1", 69_998..70_006).unwrap();
        assert_eq!(fetched, &seq[69_998..70_006]);
    }
}
