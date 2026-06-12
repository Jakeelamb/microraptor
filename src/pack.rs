//! Base and quality packing utilities.
//!
//! Bases are packed four per byte using A=0, C=1, G=2, and T=3. Ambiguous or
//! non-canonical bases are represented by zero bits in the packed byte stream
//! and a set bit in the separate ambiguity mask. Quality helpers interpret
//! input as Phred+33 FASTQ qualities.

#[cfg(all(feature = "simd", target_arch = "x86_64"))]
use std::arch::x86_64::{
    __m256i, _mm256_add_epi64, _mm256_cmpgt_epi8, _mm256_loadu_si256, _mm256_max_epu8,
    _mm256_min_epu8, _mm256_movemask_epi8, _mm256_or_si256, _mm256_sad_epu8, _mm256_set1_epi8,
    _mm256_setzero_si256, _mm256_storeu_si256, _mm256_sub_epi8,
};
use std::fmt;
use std::io::Read;
#[cfg(feature = "simd")]
use std::simd::{
    Select, Simd,
    cmp::{SimdPartialEq, SimdPartialOrd},
};

use crate::scan::scan_newlines;
use crate::{FastqConfig, FastqError, FastqPosition, Result as FastqResult};

/// Summary of a packed DNA sequence.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BaseSummary {
    /// Number of bases observed.
    pub len: usize,
    /// Number of canonical A bases.
    pub a: usize,
    /// Number of canonical C bases.
    pub c: usize,
    /// Number of canonical G bases.
    pub g: usize,
    /// Number of canonical T bases.
    pub t: usize,
    /// Number of ambiguous or non-canonical bases.
    pub n: usize,
}

impl BaseSummary {
    /// Return the number of C or G bases.
    pub fn gc_bases(self) -> usize {
        self.c + self.g
    }

    /// Return the number of A/C/G/T bases.
    pub fn canonical_bases(self) -> usize {
        self.a + self.c + self.g + self.t
    }

    /// Return true when no bases were observed.
    pub fn is_empty(self) -> bool {
        self.len == 0
    }
}

/// Packed sequence bytes plus one bit per ambiguous/non-ACGT base.
///
/// Bases are packed four per byte with base 0 in the least significant two
/// bits. Canonical bases use A=0, C=1, G=2, T=3. Masked bases store 0 in the
/// packed stream and set the corresponding bit in `n_mask`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackedSequence {
    /// Two-bit packed base bytes, four bases per byte.
    pub bases: Vec<u8>,
    /// Ambiguity bit mask, one bit per input base.
    pub n_mask: Vec<u8>,
    /// Base counts for the original sequence.
    pub summary: BaseSummary,
}

impl PackedSequence {
    /// Return the number of original bases.
    pub fn len(&self) -> usize {
        self.summary.len
    }

    /// Return true when the original sequence was empty.
    pub fn is_empty(&self) -> bool {
        self.summary.is_empty()
    }
}

/// Decoded base value from a packed sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackedBase {
    /// Canonical A.
    A,
    /// Canonical C.
    C,
    /// Canonical G.
    G,
    /// Canonical T.
    T,
    /// Ambiguous or non-canonical base.
    N,
}

/// Output buffer involved in a packing error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackBuffer {
    /// Two-bit base buffer.
    Bases,
    /// Ambiguity mask buffer.
    NMask,
    /// Quality-bin output buffer.
    QualityBins,
}

/// Error returned by base or quality packing helpers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackError {
    /// A caller-supplied output buffer was too small.
    OutputTooSmall {
        /// Buffer that was too small.
        buffer: PackBuffer,
        /// Required number of bytes.
        needed: usize,
        /// Provided number of bytes.
        provided: usize,
    },
    /// A quality byte was outside the printable Phred+33 range.
    InvalidQuality {
        /// Offset within the supplied quality slice.
        offset: usize,
        /// Invalid byte value.
        byte: u8,
    },
    /// Quality thresholds were not sorted in ascending order.
    UnsortedQualityThresholds {
        /// Threshold index that violates sorted order.
        index: usize,
    },
    /// More quality thresholds were provided than fit in a `u8` bin.
    TooManyQualityThresholds {
        /// Number of thresholds supplied by the caller.
        count: usize,
    },
}

/// Phred+33 quality summary.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct QualitySummary {
    /// Number of quality values observed.
    pub len: usize,
    /// Minimum Phred value, or `None` when empty.
    pub min_phred: Option<u8>,
    /// Maximum Phred value, or `None` when empty.
    pub max_phred: Option<u8>,
    /// Sum of Phred values.
    pub sum_phred: u64,
    /// Number of bases with Phred score at least 20.
    pub q20_bases: usize,
    /// Number of bases with Phred score at least 30.
    pub q30_bases: usize,
}

impl QualitySummary {
    /// Return the arithmetic mean Phred score, or `None` when empty.
    pub fn mean_phred(self) -> Option<f64> {
        if self.len == 0 {
            None
        } else {
            Some(self.sum_phred as f64 / self.len as f64)
        }
    }

    /// Return true when no qualities were observed.
    pub fn is_empty(self) -> bool {
        self.len == 0
    }

    fn observe(&mut self, phred: u8) {
        self.len += 1;
        self.min_phred = Some(self.min_phred.map_or(phred, |min| min.min(phred)));
        self.max_phred = Some(self.max_phred.map_or(phred, |max| max.max(phred)));
        self.sum_phred += u64::from(phred);
        self.q20_bases += usize::from(phred >= 20);
        self.q30_bases += usize::from(phred >= 30);
    }
}

/// Combined base and quality summary for one FASTQ record.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PackedRecordSummary {
    /// Packed base summary.
    pub bases: BaseSummary,
    /// Phred+33 quality summary.
    pub qualities: QualitySummary,
}

/// Borrowed view of one trusted packed FASTQ record.
///
/// The packed buffers are reused by streaming pack paths. Callers must consume
/// or copy them before the next callback invocation.
#[derive(Debug, Clone, Copy)]
pub struct TrustedPackedRecord<'a> {
    /// Record name without the leading `@`.
    pub name: &'a [u8],
    /// Original sequence bytes.
    pub seq: &'a [u8],
    /// Original quality bytes.
    pub qual: &'a [u8],
    /// Two-bit packed base bytes.
    pub bases: &'a [u8],
    /// Ambiguity bit mask.
    pub n_mask: &'a [u8],
    /// Base and quality summary for the record.
    pub summary: PackedRecordSummary,
}

/// Borrowed view of one trusted packed R1/R2 pair.
#[derive(Debug, Clone, Copy)]
pub struct TrustedPackedPair<'a> {
    /// First mate.
    pub first: TrustedPackedRecord<'a>,
    /// Second mate.
    pub second: TrustedPackedRecord<'a>,
}

/// Selected implementation family for base packing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackKernel {
    /// Portable scalar implementation.
    Scalar,
    /// Nightly portable-SIMD implementation.
    PortableSimd,
    /// x86-64 AVX2 implementation.
    Avx2,
}

/// Return the base-packing kernel selected for this build and host.
pub fn selected_pack_kernel() -> PackKernel {
    select_pack_kernel()
}

#[cfg(feature = "simd")]
fn select_pack_kernel() -> PackKernel {
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        return PackKernel::Avx2;
    }
    PackKernel::PortableSimd
}

#[cfg(not(feature = "simd"))]
fn select_pack_kernel() -> PackKernel {
    PackKernel::Scalar
}

/// Progress notification for one trusted pack slab.
#[derive(Debug, Clone, Copy)]
pub struct TrustedPackSlab {
    /// Number of complete records emitted from the slab.
    pub records: u64,
}

/// Callback interface for trusted streaming pack paths.
pub trait TrustedPackSink {
    /// Observe one packed record.
    fn record(&mut self, record: TrustedPackedRecord<'_>) -> FastqResult<()>;

    /// Observe slab-level progress after records have been emitted.
    fn slab(&mut self, _slab: TrustedPackSlab) -> FastqResult<()> {
        Ok(())
    }
}

impl<F> TrustedPackSink for F
where
    F: FnMut(TrustedPackedRecord<'_>) -> FastqResult<()>,
{
    fn record(&mut self, record: TrustedPackedRecord<'_>) -> FastqResult<()> {
        self(record)
    }
}

impl fmt::Display for PackError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::OutputTooSmall {
                buffer,
                needed,
                provided,
            } => write!(
                f,
                "{buffer:?} output too small: need {needed} bytes, got {provided}"
            ),
            Self::InvalidQuality { offset, byte } => {
                write!(f, "invalid Phred+33 quality byte {byte} at offset {offset}")
            }
            Self::UnsortedQualityThresholds { index } => {
                write!(f, "quality threshold at index {index} is not sorted")
            }
            Self::TooManyQualityThresholds { count } => {
                write!(f, "too many quality thresholds: {count}")
            }
        }
    }
}

impl std::error::Error for PackError {}

/// Pack all complete records from an in-memory trusted four-line FASTQ buffer.
///
/// This path validates record shape and quality length, but it assumes ordinary
/// four-line FASTQ and does not support multiline sequence or quality fields.
pub fn pack_trusted_fastq(
    input: &[u8],
    on_record: impl FnMut(TrustedPackedRecord<'_>) -> FastqResult<()>,
) -> FastqResult<()> {
    pack_trusted_fastq_sink(input, on_record)
}

/// Pack all complete records from an in-memory trusted FASTQ buffer into a sink.
pub fn pack_trusted_fastq_sink(input: &[u8], mut sink: impl TrustedPackSink) -> FastqResult<()> {
    let mut bases = Vec::new();
    let mut n_mask = Vec::new();
    let mut newlines = Vec::with_capacity(input.len() / 48);
    let slab = pack_trusted_fastq_slab(
        input,
        SlabContext {
            base_offset: 0,
            first_record_index: 0,
            eof: true,
        },
        &mut newlines,
        &mut bases,
        &mut n_mask,
        &mut sink,
    )?;
    debug_assert_eq!(slab.next_start, input.len());
    Ok(())
}

/// Stream trusted four-line FASTQ from a reader and invoke a callback per record.
pub fn pack_trusted_fastq_read<R: Read>(
    mut reader: R,
    config: FastqConfig,
    on_record: impl FnMut(TrustedPackedRecord<'_>) -> FastqResult<()>,
) -> FastqResult<()> {
    pack_trusted_fastq_read_sink(&mut reader, config, on_record)
}

/// Stream trusted four-line FASTQ from a reader into a sink.
pub fn pack_trusted_fastq_read_sink<R: Read>(
    mut reader: R,
    config: FastqConfig,
    mut sink: impl TrustedPackSink,
) -> FastqResult<()> {
    pack_trusted_fastq_read_sink_with_kernel(
        &mut reader,
        config,
        &mut sink,
        TrustedScanKernel::Offset,
    )
}

/// Pack an in-memory trusted FASTQ buffer with the direct scanner.
///
/// This is a lower-level alternative used for benchmark and implementation
/// comparison against the default newline-offset scanner.
pub fn pack_trusted_fastq_direct(
    input: &[u8],
    on_record: impl FnMut(TrustedPackedRecord<'_>) -> FastqResult<()>,
) -> FastqResult<()> {
    pack_trusted_fastq_direct_sink(input, on_record)
}

/// Pack an in-memory trusted FASTQ buffer with the direct scanner into a sink.
pub fn pack_trusted_fastq_direct_sink(
    input: &[u8],
    mut sink: impl TrustedPackSink,
) -> FastqResult<()> {
    let mut bases = Vec::new();
    let mut n_mask = Vec::new();
    let slab = pack_trusted_fastq_direct_slab(
        input,
        SlabContext {
            base_offset: 0,
            first_record_index: 0,
            eof: true,
        },
        &mut bases,
        &mut n_mask,
        &mut sink,
    )?;
    debug_assert_eq!(slab.next_start, input.len());
    Ok(())
}

/// Stream trusted FASTQ from a reader using the direct scanner.
pub fn pack_trusted_fastq_read_direct<R: Read>(
    mut reader: R,
    config: FastqConfig,
    on_record: impl FnMut(TrustedPackedRecord<'_>) -> FastqResult<()>,
) -> FastqResult<()> {
    pack_trusted_fastq_read_direct_sink(&mut reader, config, on_record)
}

/// Stream trusted FASTQ from a reader using the direct scanner into a sink.
pub fn pack_trusted_fastq_read_direct_sink<R: Read>(
    mut reader: R,
    config: FastqConfig,
    mut sink: impl TrustedPackSink,
) -> FastqResult<()> {
    pack_trusted_fastq_read_sink_with_kernel(
        &mut reader,
        config,
        &mut sink,
        TrustedScanKernel::Direct,
    )
}

fn pack_trusted_fastq_read_sink_with_kernel<R: Read>(
    reader: R,
    config: FastqConfig,
    sink: &mut impl TrustedPackSink,
    kernel: TrustedScanKernel,
) -> FastqResult<()> {
    let mut reader = TrustedFastqPackReader::new(reader, config, kernel);
    while reader.next_slab(sink)? {}
    Ok(())
}

struct TrustedFastqPackReader<R> {
    reader: R,
    slab_size: usize,
    buf: Vec<u8>,
    len: usize,
    eof: bool,
    base_offset: u64,
    record_index: u64,
    newlines: Vec<usize>,
    bases: Vec<u8>,
    n_mask: Vec<u8>,
    kernel: TrustedScanKernel,
}

impl<R: Read> TrustedFastqPackReader<R> {
    fn new(reader: R, config: FastqConfig, kernel: TrustedScanKernel) -> Self {
        let slab_size = config.slab_size.max(1024);
        Self {
            reader,
            slab_size,
            buf: vec![0_u8; slab_size],
            len: 0,
            eof: false,
            base_offset: 0,
            record_index: 0,
            newlines: Vec::with_capacity(slab_size / 48),
            bases: Vec::new(),
            n_mask: Vec::new(),
            kernel,
        }
    }

    fn next_slab(&mut self, sink: &mut impl TrustedPackSink) -> FastqResult<bool> {
        if self.eof && self.len == 0 {
            return Ok(false);
        }

        while !self.eof && self.len < self.slab_size {
            let n = self.reader.read(&mut self.buf[self.len..self.slab_size])?;
            if n == 0 {
                self.eof = true;
                break;
            }
            self.len += n;
        }

        let context = SlabContext {
            base_offset: self.base_offset,
            first_record_index: self.record_index,
            eof: self.eof,
        };
        let slab = match self.kernel {
            TrustedScanKernel::Offset => pack_trusted_fastq_slab(
                &self.buf[..self.len],
                context,
                &mut self.newlines,
                &mut self.bases,
                &mut self.n_mask,
                sink,
            )?,
            TrustedScanKernel::Direct => pack_trusted_fastq_direct_slab(
                &self.buf[..self.len],
                context,
                &mut self.bases,
                &mut self.n_mask,
                sink,
            )?,
        };
        if slab.records != 0 {
            sink.slab(TrustedPackSlab {
                records: slab.records,
            })?;
        }
        self.record_index += slab.records;

        if slab.next_start == self.len {
            self.base_offset += self.len as u64;
            self.len = 0;
        } else {
            let carry = self.len - slab.next_start;
            if slab.next_start == 0 && carry == self.slab_size && !self.eof {
                return Err(FastqError::RecordTooLarge {
                    slab_size: self.slab_size,
                });
            }
            self.buf.copy_within(slab.next_start..self.len, 0);
            self.base_offset += slab.next_start as u64;
            self.len = carry;
        }

        if self.eof {
            if self.len == 0 {
                return Ok(slab.records != 0);
            }
            return Err(FastqError::RecordTooLarge {
                slab_size: self.slab_size,
            });
        }

        Ok(true)
    }
}

#[derive(Clone, Copy)]
enum TrustedScanKernel {
    Offset,
    Direct,
}

/// Stream ordered R1/R2 FASTQ inputs and emit trusted packed pairs.
///
/// Pair validation uses the supplied [`crate::PairValidation`] mode. This path
/// does not synchronize reordered mates; it expects the two streams to be in
/// lockstep order.
pub fn pack_trusted_paired_fastq_read<R1: Read, R2: Read>(
    first: R1,
    second: R2,
    config: FastqConfig,
    pair_validation: crate::PairValidation,
    mut on_pair: impl FnMut(TrustedPackedPair<'_>) -> FastqResult<()>,
) -> FastqResult<()> {
    let mut first_reader = TrustedFastqLineReader::new(first, config.clone());
    let mut second_reader = TrustedFastqLineReader::new(second, config);
    let mut first_bases = Vec::new();
    let mut first_n_mask = Vec::new();
    let mut second_bases = Vec::new();
    let mut second_n_mask = Vec::new();
    let mut pair_index = 0_u64;

    loop {
        let first_count = first_reader.available_records()?;
        let second_count = second_reader.available_records()?;

        if first_count == 0 || second_count == 0 {
            if first_count == 0
                && second_count == 0
                && first_reader.is_done()
                && second_reader.is_done()
            {
                return Ok(());
            }
            if (first_count == 0 && first_reader.is_done())
                || (second_count == 0 && second_reader.is_done())
            {
                return Err(FastqError::Format(
                    "paired FASTQ inputs have different record counts".into(),
                ));
            }
            continue;
        }

        let pairs = first_count.min(second_count);
        for index in 0..pairs {
            let first = first_reader.record_lines(index);
            let second = second_reader.record_lines(index);
            let first_record_index = first_reader.record_index + index as u64;
            let second_record_index = second_reader.record_index + index as u64;

            validate_trusted_record_lines(
                first.name,
                first.seq,
                first.plus,
                first.qual,
                first_reader.base_offset,
                first_record_index,
            )?;
            validate_trusted_record_lines(
                second.name,
                second.seq,
                second.plus,
                second.qual,
                second_reader.base_offset,
                second_record_index,
            )?;

            if pair_validation != crate::PairValidation::None
                && !trusted_pair_ids_match(first.name.bytes, second.name.bytes, pair_validation)
            {
                return Err(FastqError::FormatAt {
                    message: "paired FASTQ record identifiers do not match".into(),
                    position: FastqPosition::new(0, pair_index, 0),
                });
            }

            let first_summary = pack_bases_and_summarize_qualities_into(
                first.seq.bytes,
                first.qual.bytes,
                &mut first_bases,
                &mut first_n_mask,
            )
            .map_err(|err| {
                format_fastq_at(
                    err.to_string(),
                    first_reader.base_offset,
                    first.qual.start,
                    first_record_index,
                    3,
                )
            })?;
            let second_summary = pack_bases_and_summarize_qualities_into(
                second.seq.bytes,
                second.qual.bytes,
                &mut second_bases,
                &mut second_n_mask,
            )
            .map_err(|err| {
                format_fastq_at(
                    err.to_string(),
                    second_reader.base_offset,
                    second.qual.start,
                    second_record_index,
                    3,
                )
            })?;

            on_pair(TrustedPackedPair {
                first: TrustedPackedRecord {
                    name: first.name.bytes,
                    seq: first.seq.bytes,
                    qual: first.qual.bytes,
                    bases: &first_bases,
                    n_mask: &first_n_mask,
                    summary: first_summary,
                },
                second: TrustedPackedRecord {
                    name: second.name.bytes,
                    seq: second.seq.bytes,
                    qual: second.qual.bytes,
                    bases: &second_bases,
                    n_mask: &second_n_mask,
                    summary: second_summary,
                },
            })?;
            pair_index += 1;
        }

        first_reader.consume_records(pairs)?;
        second_reader.consume_records(pairs)?;
    }
}

struct TrustedFastqLineReader<R> {
    reader: R,
    slab_size: usize,
    buf: Vec<u8>,
    len: usize,
    eof: bool,
    base_offset: u64,
    record_index: u64,
    newlines: Vec<usize>,
    line_count: usize,
    available_records: usize,
    scanned: bool,
}

impl<R: Read> TrustedFastqLineReader<R> {
    fn new(reader: R, config: FastqConfig) -> Self {
        let slab_size = config.slab_size.max(1024);
        Self {
            reader,
            slab_size,
            buf: vec![0_u8; slab_size],
            len: 0,
            eof: false,
            base_offset: 0,
            record_index: 0,
            newlines: Vec::with_capacity(slab_size / 48),
            line_count: 0,
            available_records: 0,
            scanned: false,
        }
    }

    fn available_records(&mut self) -> FastqResult<usize> {
        if self.scanned {
            return Ok(self.available_records);
        }
        if self.eof && self.len == 0 {
            self.available_records = 0;
            self.line_count = 0;
            self.scanned = true;
            return Ok(0);
        }

        while !self.eof && self.len < self.slab_size {
            let n = self.reader.read(&mut self.buf[self.len..self.slab_size])?;
            if n == 0 {
                self.eof = true;
                break;
            }
            self.len += n;
        }

        scan_newlines(&self.buf[..self.len], &mut self.newlines);
        let has_final_line = self.eof
            && self
                .newlines
                .last()
                .map_or(self.len != 0, |&nl| nl + 1 < self.len);
        self.line_count = self.newlines.len() + usize::from(has_final_line);
        if self.eof && !self.line_count.is_multiple_of(4) {
            let record_index = self.record_index + (self.line_count / 4) as u64;
            return Err(format_fastq_at(
                "truncated FASTQ record",
                self.base_offset,
                line_start(&self.newlines, (self.line_count / 4) * 4),
                record_index,
                (self.line_count % 4) as u8,
            ));
        }

        let complete_lines = (self.line_count / 4) * 4;
        self.available_records = complete_lines / 4;
        self.scanned = true;
        if self.available_records == 0 && self.len == self.slab_size && !self.eof {
            return Err(FastqError::RecordTooLarge {
                slab_size: self.slab_size,
            });
        }
        Ok(self.available_records)
    }

    fn is_done(&self) -> bool {
        self.eof && self.len == 0
    }

    fn record_lines(&self, index: usize) -> TrustedRecordLines<'_> {
        let line = index * 4;
        TrustedRecordLines {
            name: line_range(&self.buf[..self.len], &self.newlines, line),
            seq: line_range(&self.buf[..self.len], &self.newlines, line + 1),
            plus: line_range(&self.buf[..self.len], &self.newlines, line + 2),
            qual: line_range(&self.buf[..self.len], &self.newlines, line + 3),
        }
    }

    fn consume_records(&mut self, records: usize) -> FastqResult<()> {
        if records == 0 {
            return Ok(());
        }
        let consumed_lines = records * 4;
        let next_start = if self.eof && consumed_lines == self.line_count {
            self.len
        } else {
            line_start(&self.newlines, consumed_lines)
        };
        if next_start == self.len {
            self.base_offset += self.len as u64;
            self.len = 0;
        } else {
            let carry = self.len - next_start;
            if next_start == 0 && carry == self.slab_size && !self.eof {
                return Err(FastqError::RecordTooLarge {
                    slab_size: self.slab_size,
                });
            }
            self.buf.copy_within(next_start..self.len, 0);
            self.base_offset += next_start as u64;
            self.len = carry;
        }
        self.record_index += records as u64;
        self.scanned = false;
        self.line_count = 0;
        self.available_records = 0;
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct TrustedRecordLines<'a> {
    name: Line<'a>,
    seq: Line<'a>,
    plus: Line<'a>,
    qual: Line<'a>,
}

/// Return the number of bytes needed to store `base_count` two-bit bases.
pub const fn packed_base_len(base_count: usize) -> usize {
    base_count / 4 + if base_count.is_multiple_of(4) { 0 } else { 1 }
}

/// Return the number of bytes needed for a one-bit-per-item mask.
pub const fn bit_mask_len(bit_count: usize) -> usize {
    bit_count / 8 + if bit_count.is_multiple_of(8) { 0 } else { 1 }
}

const BASE_N: u8 = 4;
const BASE_LUT: [u8; 256] = base_lut();
const BASE_QUAD_STATES: usize = 5 * 5 * 5 * 5;
const BASE_QUAD_LUT: [u32; BASE_QUAD_STATES] = base_quad_lut();

const fn base_lut() -> [u8; 256] {
    let mut table = [BASE_N; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
}

const fn base_quad_lut() -> [u32; BASE_QUAD_STATES] {
    let mut table = [0_u32; BASE_QUAD_STATES];
    let mut c0 = 0_u8;
    while c0 <= BASE_N {
        let mut c1 = 0_u8;
        while c1 <= BASE_N {
            let mut c2 = 0_u8;
            while c2 <= BASE_N {
                let mut c3 = 0_u8;
                while c3 <= BASE_N {
                    let key = quad_key(c0, c1, c2, c3);
                    table[key] = quad_entry(c0, c1, c2, c3);
                    c3 += 1;
                }
                c2 += 1;
            }
            c1 += 1;
        }
        c0 += 1;
    }
    table
}

const fn quad_key(c0: u8, c1: u8, c2: u8, c3: u8) -> usize {
    c0 as usize + (c1 as usize * 5) + (c2 as usize * 25) + (c3 as usize * 125)
}

const fn quad_entry(c0: u8, c1: u8, c2: u8, c3: u8) -> u32 {
    let codes = [c0, c1, c2, c3];
    let mut packed = 0_u32;
    let mut mask = 0_u32;
    let mut counts = [0_u32; 5];
    let mut i = 0;
    while i < 4 {
        let code = codes[i];
        if code < BASE_N {
            packed |= (code as u32) << (i * 2);
            counts[code as usize] += 1;
        } else {
            mask |= 1 << i;
            counts[BASE_N as usize] += 1;
        }
        i += 1;
    }

    packed
        | (mask << 8)
        | (counts[0] << 12)
        | (counts[1] << 15)
        | (counts[2] << 18)
        | (counts[3] << 21)
        | (counts[4] << 24)
}

/// Pack a sequence into newly allocated packed-base and ambiguity buffers.
pub fn pack_bases(seq: &[u8]) -> PackedSequence {
    let mut bases = vec![0; packed_base_len(seq.len())];
    let mut n_mask = vec![0; bit_mask_len(seq.len())];
    let summary = pack_bases_exact(seq, &mut bases, &mut n_mask);
    PackedSequence {
        bases,
        n_mask,
        summary,
    }
}

/// Pack a sequence into reusable `Vec` buffers and return base counts.
pub fn pack_bases_into(seq: &[u8], bases: &mut Vec<u8>, n_mask: &mut Vec<u8>) -> BaseSummary {
    bases.clear();
    n_mask.clear();
    bases.resize(packed_base_len(seq.len()), 0);
    n_mask.resize(bit_mask_len(seq.len()), 0);
    pack_bases_exact(seq, bases, n_mask)
}

/// Pack bases and summarize Phred+33 qualities into reusable buffers.
///
/// When sequence and quality lengths match, the implementation may fuse base
/// packing and quality summary work in one pass.
#[inline]
pub fn pack_bases_and_summarize_qualities_into(
    seq: &[u8],
    qualities: &[u8],
    bases: &mut Vec<u8>,
    n_mask: &mut Vec<u8>,
) -> Result<PackedRecordSummary, PackError> {
    bases.clear();
    n_mask.clear();
    bases.resize(packed_base_len(seq.len()), 0);
    n_mask.resize(bit_mask_len(seq.len()), 0);

    if seq.len() == qualities.len() {
        pack_bases_and_qualities_exact(seq, qualities, bases, n_mask)
    } else {
        Ok(PackedRecordSummary {
            bases: pack_bases_exact(seq, bases, n_mask),
            qualities: summarize_qualities(qualities)?,
        })
    }
}

/// Pack a sequence into caller-provided fixed-size slices.
///
/// The slices may be larger than needed; only the required prefix is written.
pub fn pack_bases_into_slices(
    seq: &[u8],
    bases: &mut [u8],
    n_mask: &mut [u8],
) -> Result<BaseSummary, PackError> {
    let bases_needed = packed_base_len(seq.len());
    let mask_needed = bit_mask_len(seq.len());
    if bases.len() < bases_needed {
        return Err(PackError::OutputTooSmall {
            buffer: PackBuffer::Bases,
            needed: bases_needed,
            provided: bases.len(),
        });
    }
    if n_mask.len() < mask_needed {
        return Err(PackError::OutputTooSmall {
            buffer: PackBuffer::NMask,
            needed: mask_needed,
            provided: n_mask.len(),
        });
    }

    Ok(pack_bases_exact(
        seq,
        &mut bases[..bases_needed],
        &mut n_mask[..mask_needed],
    ))
}

/// Decode one base at `index` from packed-base and ambiguity buffers.
pub fn packed_base_at(bases: &[u8], n_mask: &[u8], index: usize) -> Option<PackedBase> {
    let base_byte = *bases.get(index / 4)?;
    if is_masked(n_mask, index)? {
        return Some(PackedBase::N);
    }

    match (base_byte >> ((index % 4) * 2)) & 0b11 {
        0 => Some(PackedBase::A),
        1 => Some(PackedBase::C),
        2 => Some(PackedBase::G),
        3 => Some(PackedBase::T),
        _ => None,
    }
}

/// Return whether `index` is marked ambiguous in an ambiguity mask.
pub fn is_masked(n_mask: &[u8], index: usize) -> Option<bool> {
    let mask_byte = *n_mask.get(index / 8)?;
    Some(((mask_byte >> (index % 8)) & 1) != 0)
}

/// Summarize a slice of Phred+33 quality bytes.
pub fn summarize_qualities(qualities: &[u8]) -> Result<QualitySummary, PackError> {
    #[cfg(all(feature = "simd", target_arch = "x86_64"))]
    if qualities.len() >= 32 && std::is_x86_feature_detected!("avx2") {
        return unsafe { summarize_qualities_avx2(qualities) };
    }
    #[cfg(feature = "simd")]
    if qualities.len() >= 32 {
        return summarize_qualities_simd(qualities);
    }
    let mut summary = QualityAccumulator::default();
    for (offset, &byte) in qualities.iter().enumerate() {
        summary.observe(byte, offset)?;
    }
    Ok(summary.finish())
}

#[cfg(all(feature = "simd", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn summarize_qualities_avx2(qualities: &[u8]) -> Result<QualitySummary, PackError> {
    let low = _mm256_set1_epi8(33);
    let high = _mm256_set1_epi8(126);
    let offset = _mm256_set1_epi8(33);
    let q20 = _mm256_set1_epi8(52);
    let q30 = _mm256_set1_epi8(62);
    let zero = _mm256_setzero_si256();
    let mut min_phred = _mm256_set1_epi8(93);
    let mut max_phred = _mm256_setzero_si256();
    let mut sum_phred = _mm256_setzero_si256();

    let mut summary = QualityAccumulator::default();
    let mut i = 0;
    while i + 32 <= qualities.len() {
        let bytes = unsafe { _mm256_loadu_si256(qualities.as_ptr().add(i).cast::<__m256i>()) };
        let too_low = _mm256_cmpgt_epi8(low, bytes);
        let too_high = _mm256_cmpgt_epi8(bytes, high);
        let invalid = _mm256_movemask_epi8(_mm256_or_si256(too_low, too_high));
        if invalid != 0 {
            let offset = invalid.trailing_zeros() as usize;
            return Err(PackError::InvalidQuality {
                offset: i + offset,
                byte: qualities[i + offset],
            });
        }

        summary.q20_bases +=
            _mm256_movemask_epi8(_mm256_cmpgt_epi8(bytes, q20)).count_ones() as usize;
        summary.q30_bases +=
            _mm256_movemask_epi8(_mm256_cmpgt_epi8(bytes, q30)).count_ones() as usize;
        let phreds = _mm256_sub_epi8(bytes, offset);
        min_phred = _mm256_min_epu8(min_phred, phreds);
        max_phred = _mm256_max_epu8(max_phred, phreds);
        sum_phred = _mm256_add_epi64(sum_phred, _mm256_sad_epu8(phreds, zero));
        summary.len += 32;
        i += 32;
    }

    unsafe { finish_avx2_quality_vectors(&mut summary, min_phred, max_phred, sum_phred) };

    while i < qualities.len() {
        summary.observe(qualities[i], i)?;
        i += 1;
    }

    Ok(summary.finish())
}

#[cfg(all(feature = "simd", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn finish_avx2_quality_vectors(
    summary: &mut QualityAccumulator,
    min_phred: __m256i,
    max_phred: __m256i,
    sum_phred: __m256i,
) {
    if summary.len == 0 {
        return;
    }

    let mut min_lanes = [0_u8; 32];
    let mut max_lanes = [0_u8; 32];
    let mut sum_lanes = [0_u64; 4];
    unsafe {
        _mm256_storeu_si256(min_lanes.as_mut_ptr().cast::<__m256i>(), min_phred);
        _mm256_storeu_si256(max_lanes.as_mut_ptr().cast::<__m256i>(), max_phred);
        _mm256_storeu_si256(sum_lanes.as_mut_ptr().cast::<__m256i>(), sum_phred);
    }
    for &phred in &min_lanes {
        summary.min_phred = summary.min_phred.min(phred);
    }
    for &phred in &max_lanes {
        summary.max_phred = summary.max_phred.max(phred);
    }
    summary.sum_phred = sum_lanes.iter().copied().sum();
}

#[cfg(feature = "simd")]
fn summarize_qualities_simd(qualities: &[u8]) -> Result<QualitySummary, PackError> {
    const LANES: usize = 32;
    type Chunk = Simd<u8, LANES>;

    let low = Chunk::splat(33);
    let high = Chunk::splat(126);
    let offset = Chunk::splat(33);
    let q20 = Chunk::splat(20);
    let q30 = Chunk::splat(30);

    let mut summary = QualityAccumulator::default();
    let mut i = 0;
    while i + LANES <= qualities.len() {
        let bytes = Chunk::from_slice(&qualities[i..i + LANES]);
        if (bytes.simd_lt(low) | bytes.simd_gt(high)).any() {
            let mut j = 0;
            while j < LANES {
                phred33(qualities[i + j], i + j)?;
                j += 1;
            }
        }
        let phreds = bytes - offset;
        for phred in phreds.to_array() {
            summary.min_phred = summary.min_phred.min(phred);
            summary.max_phred = summary.max_phred.max(phred);
            summary.sum_phred += u64::from(phred);
        }
        summary.len += LANES;
        summary.q20_bases += phreds.simd_ge(q20).to_bitmask().count_ones() as usize;
        summary.q30_bases += phreds.simd_ge(q30).to_bitmask().count_ones() as usize;
        i += LANES;
    }

    while i < qualities.len() {
        summary.observe(qualities[i], i)?;
        i += 1;
    }

    Ok(summary.finish())
}

/// Bin Phred+33 qualities into threshold indexes.
///
/// Thresholds are lower bounds for the next bin. For thresholds `[10, 20, 30]`,
/// Phred 0-9 maps to 0, 10-19 maps to 1, 20-29 maps to 2, and 30+ maps to 3.
pub fn bin_qualities_into(
    qualities: &[u8],
    thresholds: &[u8],
    out: &mut Vec<u8>,
) -> Result<QualitySummary, PackError> {
    validate_thresholds(thresholds)?;
    out.clear();
    out.reserve(qualities.len());

    let mut summary = QualitySummary::default();
    for (offset, &byte) in qualities.iter().enumerate() {
        let phred = phred33(byte, offset)?;
        summary.observe(phred);
        out.push(quality_bin(phred, thresholds));
    }
    Ok(summary)
}

/// Bin Phred+33 qualities into a caller-provided output slice.
///
/// See [`bin_qualities_into`] for threshold semantics.
pub fn bin_qualities_into_slice(
    qualities: &[u8],
    thresholds: &[u8],
    out: &mut [u8],
) -> Result<QualitySummary, PackError> {
    validate_thresholds(thresholds)?;
    if out.len() < qualities.len() {
        return Err(PackError::OutputTooSmall {
            buffer: PackBuffer::QualityBins,
            needed: qualities.len(),
            provided: out.len(),
        });
    }

    let mut summary = QualitySummary::default();
    for (offset, &byte) in qualities.iter().enumerate() {
        let phred = phred33(byte, offset)?;
        summary.observe(phred);
        out[offset] = quality_bin(phred, thresholds);
    }
    Ok(summary)
}

fn pack_trusted_fastq_slab(
    input: &[u8],
    context: SlabContext,
    newlines: &mut Vec<usize>,
    bases: &mut Vec<u8>,
    n_mask: &mut Vec<u8>,
    sink: &mut impl TrustedPackSink,
) -> FastqResult<SlabResult> {
    scan_newlines(input, newlines);
    let mut records = 0_u64;

    let has_final_line = context.eof
        && newlines
            .last()
            .map_or(!input.is_empty(), |&nl| nl + 1 < input.len());
    let line_count = newlines.len() + usize::from(has_final_line);
    if context.eof && !line_count.is_multiple_of(4) {
        let record_index = context.first_record_index + (line_count / 4) as u64;
        return Err(format_fastq_at(
            "truncated FASTQ record",
            context.base_offset,
            line_start(newlines, (line_count / 4) * 4),
            record_index,
            (line_count % 4) as u8,
        ));
    }
    let complete_lines = (line_count / 4) * 4;

    for line in (0..complete_lines).step_by(4) {
        let record_index = context.first_record_index + (line / 4) as u64;
        let name = line_range(input, newlines, line);
        let seq = line_range(input, newlines, line + 1);
        let plus = line_range(input, newlines, line + 2);
        let qual = line_range(input, newlines, line + 3);

        observe_trusted_packed_record(
            name,
            seq,
            plus,
            qual,
            context.base_offset,
            record_index,
            bases,
            n_mask,
            sink,
        )?;
        records += 1;
    }

    if complete_lines == line_count && context.eof {
        Ok(SlabResult {
            next_start: input.len(),
            records,
        })
    } else {
        let next_start = line_start(newlines, complete_lines);
        if complete_lines == line_count && next_start == input.len() {
            return Ok(SlabResult {
                next_start: input.len(),
                records,
            });
        }
        Ok(SlabResult {
            next_start,
            records,
        })
    }
}

fn pack_trusted_fastq_direct_slab(
    input: &[u8],
    context: SlabContext,
    bases: &mut Vec<u8>,
    n_mask: &mut Vec<u8>,
    sink: &mut impl TrustedPackSink,
) -> FastqResult<SlabResult> {
    let mut cursor = 0;
    let mut records = 0_u64;

    while cursor < input.len() {
        let record_start = cursor;
        let Some(name) = direct_line(input, &mut cursor, context.eof) else {
            return Ok(SlabResult {
                next_start: record_start,
                records,
            });
        };
        let Some(seq) = direct_line(input, &mut cursor, context.eof) else {
            return incomplete_or_truncated_direct(input, context, record_start, records, 1);
        };
        let Some(plus) = direct_line(input, &mut cursor, context.eof) else {
            return incomplete_or_truncated_direct(input, context, record_start, records, 2);
        };
        let Some(qual) = direct_line(input, &mut cursor, context.eof) else {
            return incomplete_or_truncated_direct(input, context, record_start, records, 3);
        };

        observe_trusted_packed_record(
            name,
            seq,
            plus,
            qual,
            context.base_offset,
            context.first_record_index + records,
            bases,
            n_mask,
            sink,
        )?;
        records += 1;
    }

    Ok(SlabResult {
        next_start: input.len(),
        records,
    })
}

fn incomplete_or_truncated_direct(
    input: &[u8],
    context: SlabContext,
    record_start: usize,
    records: u64,
    line_index: u8,
) -> FastqResult<SlabResult> {
    if context.eof {
        Err(format_fastq_at(
            "truncated FASTQ record",
            context.base_offset,
            input.len(),
            context.first_record_index + records,
            line_index,
        ))
    } else {
        Ok(SlabResult {
            next_start: record_start,
            records,
        })
    }
}

fn direct_line<'a>(input: &'a [u8], cursor: &mut usize, eof: bool) -> Option<Line<'a>> {
    let start = *cursor;
    if start >= input.len() {
        return None;
    }

    let mut end = start;
    while end < input.len() && input[end] != b'\n' {
        end += 1;
    }
    if end == input.len() && !eof {
        return None;
    }

    *cursor = if end < input.len() { end + 1 } else { end };
    let end = trim_cr_end(input, start, end);
    Some(Line {
        bytes: &input[start..end],
        start,
    })
}

#[derive(Clone, Copy)]
struct SlabContext {
    base_offset: u64,
    first_record_index: u64,
    eof: bool,
}

#[derive(Clone, Copy)]
struct SlabResult {
    next_start: usize,
    records: u64,
}

#[allow(clippy::too_many_arguments)]
fn observe_trusted_packed_record(
    name: Line<'_>,
    seq: Line<'_>,
    plus: Line<'_>,
    qual: Line<'_>,
    base_offset: u64,
    record_index: u64,
    bases: &mut Vec<u8>,
    n_mask: &mut Vec<u8>,
    sink: &mut impl TrustedPackSink,
) -> FastqResult<()> {
    validate_trusted_record_lines(name, seq, plus, qual, base_offset, record_index)?;

    let summary = pack_bases_and_summarize_qualities_into(seq.bytes, qual.bytes, bases, n_mask)
        .map_err(|err| {
            format_fastq_at(err.to_string(), base_offset, qual.start, record_index, 3)
        })?;
    sink.record(TrustedPackedRecord {
        name: name.bytes,
        seq: seq.bytes,
        qual: qual.bytes,
        bases: &bases[..],
        n_mask: &n_mask[..],
        summary,
    })
}

fn validate_trusted_record_lines(
    name: Line<'_>,
    seq: Line<'_>,
    plus: Line<'_>,
    qual: Line<'_>,
    base_offset: u64,
    record_index: u64,
) -> FastqResult<()> {
    if name.bytes.first() != Some(&b'@') {
        return Err(format_fastq_at(
            "header must start with `@`",
            base_offset,
            name.start,
            record_index,
            0,
        ));
    }
    if plus.bytes.first() != Some(&b'+') {
        return Err(format_fastq_at(
            "plus line must start with `+`",
            base_offset,
            plus.start,
            record_index,
            2,
        ));
    }
    if seq.len() != qual.len() {
        return Err(format_fastq_at(
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
    Ok(())
}

fn trusted_pair_ids_match(
    first_name: &[u8],
    second_name: &[u8],
    mode: crate::PairValidation,
) -> bool {
    match mode {
        crate::PairValidation::None => true,
        crate::PairValidation::FastSlash => fast_slash_pair_ids_match(first_name, second_name)
            .unwrap_or_else(|| normalized_pair_id(first_name) == normalized_pair_id(second_name)),
        crate::PairValidation::Full => {
            normalized_pair_id(first_name) == normalized_pair_id(second_name)
        }
    }
}

fn fast_slash_pair_ids_match(first_name: &[u8], second_name: &[u8]) -> Option<bool> {
    let first = first_name.strip_prefix(b"@").unwrap_or(first_name);
    let second = second_name.strip_prefix(b"@").unwrap_or(second_name);
    let first_end = token_end(first);
    let second_end = token_end(second);
    let first = &first[..first_end];
    let second = &second[..second_end];
    if first.len() < 3 || second.len() < 3 || first.len() != second.len() {
        return None;
    }
    if !first.ends_with(b"/1") || !second.ends_with(b"/2") {
        return None;
    }
    Some(first[..first.len() - 2] == second[..second.len() - 2])
}

fn normalized_pair_id(name: &[u8]) -> &[u8] {
    let name = name.strip_prefix(b"@").unwrap_or(name);
    let token = &name[..token_end(name)];
    if token.len() >= 2 && (token.ends_with(b"/1") || token.ends_with(b"/2")) {
        &token[..token.len() - 2]
    } else {
        token
    }
}

fn token_end(bytes: &[u8]) -> usize {
    let mut end = 0;
    while end < bytes.len() && !bytes[end].is_ascii_whitespace() {
        end += 1;
    }
    end
}

fn pack_bases_exact(seq: &[u8], bases: &mut [u8], n_mask: &mut [u8]) -> BaseSummary {
    let bases_needed = packed_base_len(seq.len());
    if bases_needed == 0 {
        return BaseSummary::default();
    }

    bases[..bases_needed].fill(0);
    let mask_needed = bit_mask_len(seq.len());
    if mask_needed != 0 {
        n_mask[..mask_needed].fill(0);
    }

    let mut summary = BaseSummary {
        len: seq.len(),
        ..BaseSummary::default()
    };

    let full_chunks = seq.len() / 4;
    let mut chunk_index = 0;
    let mut base_index = 0;
    #[cfg(feature = "simd")]
    while chunk_index + 4 <= full_chunks {
        let codes = base_codes_16(&seq[base_index..base_index + 16]);
        pack_code_quads(
            &codes,
            base_index,
            &mut bases[chunk_index..chunk_index + 4],
            &mut summary,
            n_mask,
        );
        chunk_index += 4;
        base_index += 16;
    }
    while chunk_index < full_chunks {
        let c0 = BASE_LUT[usize::from(seq[base_index])];
        let c1 = BASE_LUT[usize::from(seq[base_index + 1])];
        let c2 = BASE_LUT[usize::from(seq[base_index + 2])];
        let c3 = BASE_LUT[usize::from(seq[base_index + 3])];
        pack_quad_from_codes(
            c0,
            c1,
            c2,
            c3,
            base_index,
            &mut bases[chunk_index],
            &mut summary,
            n_mask,
        );
        chunk_index += 1;
        base_index += 4;
    }

    let tail_start = full_chunks * 4;
    let mut index = tail_start;
    while index < seq.len() {
        let offset = index - tail_start;
        let code = BASE_LUT[usize::from(seq[index])];
        if code < BASE_N {
            add_base_count(&mut summary, code);
            bases[full_chunks] |= code << (offset * 2);
        } else {
            summary.n += 1;
            n_mask[index / 8] |= 1 << (index % 8);
        }
        index += 1;
    }

    summary
}

fn pack_bases_and_qualities_exact(
    seq: &[u8],
    qualities: &[u8],
    bases: &mut [u8],
    n_mask: &mut [u8],
) -> Result<PackedRecordSummary, PackError> {
    let bases_needed = packed_base_len(seq.len());
    if bases_needed == 0 {
        return Ok(PackedRecordSummary::default());
    }

    bases[..bases_needed].fill(0);
    let mask_needed = bit_mask_len(seq.len());
    if mask_needed != 0 {
        n_mask[..mask_needed].fill(0);
    }

    #[cfg(all(feature = "simd", target_arch = "x86_64"))]
    if seq.len() >= 32 && std::is_x86_feature_detected!("avx2") {
        return unsafe { pack_bases_and_qualities_exact_avx2(seq, qualities, bases, n_mask) };
    }

    let mut bases_summary = BaseSummary {
        len: seq.len(),
        ..BaseSummary::default()
    };
    let mut quality_summary = QualityAccumulator::default();

    let full_chunks = seq.len() / 4;
    let mut chunk_index = 0;
    let mut base_index = 0;
    #[cfg(feature = "simd")]
    while chunk_index + 4 <= full_chunks {
        let codes = base_codes_16(&seq[base_index..base_index + 16]);
        pack_code_quads(
            &codes,
            base_index,
            &mut bases[chunk_index..chunk_index + 4],
            &mut bases_summary,
            n_mask,
        );
        let end = base_index + 16;
        while base_index < end {
            quality_summary.observe(qualities[base_index], base_index)?;
            base_index += 1;
        }
        chunk_index += 4;
    }
    while chunk_index < full_chunks {
        let c0 = BASE_LUT[usize::from(seq[base_index])];
        let c1 = BASE_LUT[usize::from(seq[base_index + 1])];
        let c2 = BASE_LUT[usize::from(seq[base_index + 2])];
        let c3 = BASE_LUT[usize::from(seq[base_index + 3])];
        pack_quad_from_codes(
            c0,
            c1,
            c2,
            c3,
            base_index,
            &mut bases[chunk_index],
            &mut bases_summary,
            n_mask,
        );
        quality_summary.observe(qualities[base_index], base_index)?;
        quality_summary.observe(qualities[base_index + 1], base_index + 1)?;
        quality_summary.observe(qualities[base_index + 2], base_index + 2)?;
        quality_summary.observe(qualities[base_index + 3], base_index + 3)?;
        chunk_index += 1;
        base_index += 4;
    }

    let tail_start = full_chunks * 4;
    let mut index = tail_start;
    while index < seq.len() {
        let offset = index - tail_start;
        let code = BASE_LUT[usize::from(seq[index])];
        if code < BASE_N {
            add_base_count(&mut bases_summary, code);
            bases[full_chunks] |= code << (offset * 2);
        } else {
            bases_summary.n += 1;
            n_mask[index / 8] |= 1 << (index % 8);
        }
        quality_summary.observe(qualities[index], index)?;
        index += 1;
    }

    Ok(PackedRecordSummary {
        bases: bases_summary,
        qualities: quality_summary.finish(),
    })
}

#[cfg(all(feature = "simd", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn pack_bases_and_qualities_exact_avx2(
    seq: &[u8],
    qualities: &[u8],
    bases: &mut [u8],
    n_mask: &mut [u8],
) -> Result<PackedRecordSummary, PackError> {
    let low = _mm256_set1_epi8(33);
    let high = _mm256_set1_epi8(126);
    let offset = _mm256_set1_epi8(33);
    let q20 = _mm256_set1_epi8(52);
    let q30 = _mm256_set1_epi8(62);
    let zero = _mm256_setzero_si256();
    let mut min_phred = _mm256_set1_epi8(93);
    let mut max_phred = _mm256_setzero_si256();
    let mut sum_phred = _mm256_setzero_si256();

    let mut bases_summary = BaseSummary {
        len: seq.len(),
        ..BaseSummary::default()
    };
    let mut quality_summary = QualityAccumulator::default();
    let mut base_index = 0;
    let mut chunk_index = 0;

    while base_index + 32 <= seq.len() {
        let block_start = base_index;
        let first_codes = base_codes_16(&seq[base_index..base_index + 16]);
        pack_code_quads(
            &first_codes,
            base_index,
            &mut bases[chunk_index..chunk_index + 4],
            &mut bases_summary,
            n_mask,
        );
        base_index += 16;
        chunk_index += 4;

        let second_codes = base_codes_16(&seq[base_index..base_index + 16]);
        pack_code_quads(
            &second_codes,
            base_index,
            &mut bases[chunk_index..chunk_index + 4],
            &mut bases_summary,
            n_mask,
        );
        base_index += 16;
        chunk_index += 4;

        let bytes =
            unsafe { _mm256_loadu_si256(qualities.as_ptr().add(block_start).cast::<__m256i>()) };
        let too_low = _mm256_cmpgt_epi8(low, bytes);
        let too_high = _mm256_cmpgt_epi8(bytes, high);
        let invalid = _mm256_movemask_epi8(_mm256_or_si256(too_low, too_high));
        if invalid != 0 {
            let offset = invalid.trailing_zeros() as usize;
            return Err(PackError::InvalidQuality {
                offset: block_start + offset,
                byte: qualities[block_start + offset],
            });
        }

        quality_summary.q20_bases +=
            _mm256_movemask_epi8(_mm256_cmpgt_epi8(bytes, q20)).count_ones() as usize;
        quality_summary.q30_bases +=
            _mm256_movemask_epi8(_mm256_cmpgt_epi8(bytes, q30)).count_ones() as usize;
        let phreds = _mm256_sub_epi8(bytes, offset);
        min_phred = _mm256_min_epu8(min_phred, phreds);
        max_phred = _mm256_max_epu8(max_phred, phreds);
        sum_phred = _mm256_add_epi64(sum_phred, _mm256_sad_epu8(phreds, zero));
        quality_summary.len += 32;
    }

    unsafe { finish_avx2_quality_vectors(&mut quality_summary, min_phred, max_phred, sum_phred) };

    let full_chunks = seq.len() / 4;
    while chunk_index < full_chunks {
        let c0 = BASE_LUT[usize::from(seq[base_index])];
        let c1 = BASE_LUT[usize::from(seq[base_index + 1])];
        let c2 = BASE_LUT[usize::from(seq[base_index + 2])];
        let c3 = BASE_LUT[usize::from(seq[base_index + 3])];
        pack_quad_from_codes(
            c0,
            c1,
            c2,
            c3,
            base_index,
            &mut bases[chunk_index],
            &mut bases_summary,
            n_mask,
        );
        quality_summary.observe(qualities[base_index], base_index)?;
        quality_summary.observe(qualities[base_index + 1], base_index + 1)?;
        quality_summary.observe(qualities[base_index + 2], base_index + 2)?;
        quality_summary.observe(qualities[base_index + 3], base_index + 3)?;
        chunk_index += 1;
        base_index += 4;
    }

    while base_index < seq.len() {
        let offset = base_index - (full_chunks * 4);
        let code = BASE_LUT[usize::from(seq[base_index])];
        if code < BASE_N {
            add_base_count(&mut bases_summary, code);
            bases[full_chunks] |= code << (offset * 2);
        } else {
            bases_summary.n += 1;
            n_mask[base_index / 8] |= 1 << (base_index % 8);
        }
        quality_summary.observe(qualities[base_index], base_index)?;
        base_index += 1;
    }

    Ok(PackedRecordSummary {
        bases: bases_summary,
        qualities: quality_summary.finish(),
    })
}

#[cfg(feature = "simd")]
fn base_codes_16(seq: &[u8]) -> [u8; 16] {
    debug_assert!(seq.len() >= 16);
    type Chunk = Simd<u8, 16>;

    let lower = Chunk::from_slice(&seq[..16]) | Chunk::splat(0x20);
    let mut codes = Chunk::splat(BASE_N);
    codes = lower
        .simd_eq(Chunk::splat(b'a'))
        .select(Chunk::splat(0), codes);
    codes = lower
        .simd_eq(Chunk::splat(b'c'))
        .select(Chunk::splat(1), codes);
    codes = lower
        .simd_eq(Chunk::splat(b'g'))
        .select(Chunk::splat(2), codes);
    codes = lower
        .simd_eq(Chunk::splat(b't'))
        .select(Chunk::splat(3), codes);
    codes.to_array()
}

#[cfg(feature = "simd")]
fn pack_code_quads(
    codes: &[u8; 16],
    base_index: usize,
    bases: &mut [u8],
    summary: &mut BaseSummary,
    n_mask: &mut [u8],
) {
    if codes_are_canonical(codes) {
        pack_canonical_code_quads(codes, bases, summary);
        return;
    }

    let mut quad = 0;
    while quad < 4 {
        let code_index = quad * 4;
        let index = base_index + code_index;
        pack_quad_from_codes(
            codes[code_index],
            codes[code_index + 1],
            codes[code_index + 2],
            codes[code_index + 3],
            index,
            &mut bases[quad],
            summary,
            n_mask,
        );
        quad += 1;
    }
}

#[cfg(feature = "simd")]
#[inline(always)]
fn codes_are_canonical(codes: &[u8; 16]) -> bool {
    let mut combined = 0_u8;
    let mut i = 0;
    while i < 16 {
        combined |= codes[i];
        i += 1;
    }
    combined < BASE_N
}

#[cfg(feature = "simd")]
#[inline(always)]
fn pack_canonical_code_quads(codes: &[u8; 16], bases: &mut [u8], summary: &mut BaseSummary) {
    let mut quad = 0;
    while quad < 4 {
        let i = quad * 4;
        bases[quad] = codes[i] | (codes[i + 1] << 2) | (codes[i + 2] << 4) | (codes[i + 3] << 6);
        quad += 1;
    }

    let mut i = 0;
    while i < 16 {
        add_base_count(summary, codes[i]);
        i += 1;
    }
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn pack_quad_from_codes(
    c0: u8,
    c1: u8,
    c2: u8,
    c3: u8,
    base_index: usize,
    base_out: &mut u8,
    summary: &mut BaseSummary,
    n_mask: &mut [u8],
) {
    apply_quad_entry(
        BASE_QUAD_LUT[quad_key(c0, c1, c2, c3)],
        base_index,
        base_out,
        summary,
        n_mask,
    );
}

#[inline(always)]
fn apply_quad_entry(
    entry: u32,
    base_index: usize,
    base_out: &mut u8,
    summary: &mut BaseSummary,
    n_mask: &mut [u8],
) {
    *base_out = entry as u8;
    let mask = ((entry >> 8) & 0x0f) as u8;
    if mask != 0 {
        n_mask[base_index / 8] |= mask << (base_index % 8);
    }
    summary.a += ((entry >> 12) & 0x07) as usize;
    summary.c += ((entry >> 15) & 0x07) as usize;
    summary.g += ((entry >> 18) & 0x07) as usize;
    summary.t += ((entry >> 21) & 0x07) as usize;
    summary.n += ((entry >> 24) & 0x07) as usize;
}

struct QualityAccumulator {
    len: usize,
    min_phred: u8,
    max_phred: u8,
    sum_phred: u64,
    q20_bases: usize,
    q30_bases: usize,
}

impl Default for QualityAccumulator {
    fn default() -> Self {
        Self {
            len: 0,
            min_phred: u8::MAX,
            max_phred: 0,
            sum_phred: 0,
            q20_bases: 0,
            q30_bases: 0,
        }
    }
}

impl QualityAccumulator {
    #[inline(always)]
    fn observe(&mut self, byte: u8, offset: usize) -> Result<(), PackError> {
        let phred = phred33(byte, offset)?;
        self.min_phred = self.min_phred.min(phred);
        self.max_phred = self.max_phred.max(phred);
        self.len += 1;
        self.sum_phred += u64::from(phred);
        self.q20_bases += usize::from(phred >= 20);
        self.q30_bases += usize::from(phred >= 30);
        Ok(())
    }

    #[inline(always)]
    fn finish(self) -> QualitySummary {
        if self.len == 0 {
            return QualitySummary::default();
        }
        QualitySummary {
            len: self.len,
            min_phred: Some(self.min_phred),
            max_phred: Some(self.max_phred),
            sum_phred: self.sum_phred,
            q20_bases: self.q20_bases,
            q30_bases: self.q30_bases,
        }
    }
}

#[inline(always)]
fn add_base_count(summary: &mut BaseSummary, code: u8) {
    match code {
        0 => summary.a += 1,
        1 => summary.c += 1,
        2 => summary.g += 1,
        3 => summary.t += 1,
        _ => unreachable!(),
    }
}

#[inline(always)]
fn phred33(byte: u8, offset: usize) -> Result<u8, PackError> {
    if (33..=126).contains(&byte) {
        Ok(byte - 33)
    } else {
        Err(PackError::InvalidQuality { offset, byte })
    }
}

fn validate_thresholds(thresholds: &[u8]) -> Result<(), PackError> {
    if thresholds.len() > usize::from(u8::MAX) {
        return Err(PackError::TooManyQualityThresholds {
            count: thresholds.len(),
        });
    }

    for (index, pair) in thresholds.windows(2).enumerate() {
        if pair[0] > pair[1] {
            return Err(PackError::UnsortedQualityThresholds { index: index + 1 });
        }
    }

    Ok(())
}

fn quality_bin(phred: u8, thresholds: &[u8]) -> u8 {
    let mut bin = 0;
    for &threshold in thresholds {
        if phred < threshold {
            break;
        }
        bin += 1;
    }
    bin
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

fn format_fastq_at(
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

    #[test]
    fn reports_required_lengths() {
        assert_eq!(packed_base_len(0), 0);
        assert_eq!(packed_base_len(1), 1);
        assert_eq!(packed_base_len(4), 1);
        assert_eq!(packed_base_len(5), 2);

        assert_eq!(bit_mask_len(0), 0);
        assert_eq!(bit_mask_len(1), 1);
        assert_eq!(bit_mask_len(8), 1);
        assert_eq!(bit_mask_len(9), 2);
    }

    #[test]
    fn packs_canonical_bases_four_per_byte() {
        let packed = pack_bases(b"ACGTacgt");
        assert_eq!(packed.bases, vec![0b1110_0100, 0b1110_0100]);
        assert_eq!(packed.n_mask, vec![0]);
        assert_eq!(
            packed.summary,
            BaseSummary {
                len: 8,
                a: 2,
                c: 2,
                g: 2,
                t: 2,
                n: 0,
            }
        );

        let decoded: Vec<_> = (0..packed.len())
            .map(|i| packed_base_at(&packed.bases, &packed.n_mask, i))
            .collect();
        assert_eq!(
            decoded,
            vec![
                Some(PackedBase::A),
                Some(PackedBase::C),
                Some(PackedBase::G),
                Some(PackedBase::T),
                Some(PackedBase::A),
                Some(PackedBase::C),
                Some(PackedBase::G),
                Some(PackedBase::T),
            ]
        );
    }

    #[test]
    fn masks_ambiguous_bases_without_rejecting_record() {
        let packed = pack_bases(b"ANGTXry");
        assert_eq!(packed.bases, vec![0b1110_0000, 0]);
        assert_eq!(packed.n_mask, vec![0b0111_0010]);
        assert_eq!(packed.summary.n, 4);
        assert_eq!(packed.summary.canonical_bases(), 3);

        assert_eq!(
            (0..packed.len())
                .map(|i| packed_base_at(&packed.bases, &packed.n_mask, i))
                .collect::<Vec<_>>(),
            vec![
                Some(PackedBase::A),
                Some(PackedBase::N),
                Some(PackedBase::G),
                Some(PackedBase::T),
                Some(PackedBase::N),
                Some(PackedBase::N),
                Some(PackedBase::N),
            ]
        );
    }

    #[test]
    fn reuses_vec_buffers() {
        let mut bases = Vec::with_capacity(16);
        let mut n_mask = Vec::with_capacity(16);
        bases.extend_from_slice(&[255; 8]);
        n_mask.extend_from_slice(&[255; 8]);

        let summary = pack_bases_into(b"AAAAA", &mut bases, &mut n_mask);
        assert_eq!(summary.len, 5);
        assert_eq!(bases, vec![0, 0]);
        assert_eq!(n_mask, vec![0]);
        assert!(bases.capacity() >= 16);
        assert!(n_mask.capacity() >= 16);
    }

    #[test]
    fn slice_pack_reports_small_buffers() {
        let mut bases = [0; 1];
        let mut n_mask = [0; 1];
        let err = pack_bases_into_slices(b"ACGTA", &mut bases, &mut n_mask).unwrap_err();
        assert_eq!(
            err,
            PackError::OutputTooSmall {
                buffer: PackBuffer::Bases,
                needed: 2,
                provided: 1,
            }
        );
    }

    #[test]
    fn summarizes_phred33_quality() {
        let summary = summarize_qualities(b"!5?I").unwrap();
        assert_eq!(
            summary,
            QualitySummary {
                len: 4,
                min_phred: Some(0),
                max_phred: Some(40),
                sum_phred: 90,
                q20_bases: 3,
                q30_bases: 2,
            }
        );
        assert_eq!(summary.mean_phred(), Some(22.5));
    }

    #[test]
    fn summarizes_long_quality_vector_reduction() {
        let qualities: Vec<u8> = (0..97).map(|i| 33 + (i % 41) as u8).collect();
        let summary = summarize_qualities(&qualities).unwrap();
        let phreds: Vec<u8> = qualities.iter().map(|&byte| byte - 33).collect();

        assert_eq!(summary.len, qualities.len());
        assert_eq!(summary.min_phred, phreds.iter().copied().min());
        assert_eq!(summary.max_phred, phreds.iter().copied().max());
        assert_eq!(
            summary.sum_phred,
            phreds.iter().map(|&phred| u64::from(phred)).sum()
        );
        assert_eq!(
            summary.q20_bases,
            phreds.iter().filter(|&&q| q >= 20).count()
        );
        assert_eq!(
            summary.q30_bases,
            phreds.iter().filter(|&&q| q >= 30).count()
        );
    }

    #[test]
    fn packs_bases_and_summarizes_qualities_together() {
        let mut bases = Vec::new();
        let mut n_mask = Vec::new();
        let summary =
            pack_bases_and_summarize_qualities_into(b"ACGTNN", b"!5?III", &mut bases, &mut n_mask)
                .unwrap();

        assert_eq!(bases, vec![0b1110_0100, 0]);
        assert_eq!(n_mask, vec![0b0011_0000]);
        assert_eq!(summary.bases.len, 6);
        assert_eq!(summary.bases.canonical_bases(), 4);
        assert_eq!(summary.bases.n, 2);
        assert_eq!(
            summary.qualities,
            QualitySummary {
                len: 6,
                min_phred: Some(0),
                max_phred: Some(40),
                sum_phred: 170,
                q20_bases: 5,
                q30_bases: 4,
            }
        );
    }

    #[test]
    fn fused_pack_handles_long_canonical_mixed_case_read() {
        let seq = b"ACGTacgtACGTacgtACGTacgtACGTacgtACGTacgt";
        let qual = vec![b'I'; seq.len()];
        let mut bases = Vec::new();
        let mut n_mask = Vec::new();
        let summary =
            pack_bases_and_summarize_qualities_into(seq, &qual, &mut bases, &mut n_mask).unwrap();

        assert_eq!(summary.bases.len, seq.len());
        assert_eq!(summary.bases.n, 0);
        assert_eq!(summary.bases.a, 10);
        assert_eq!(summary.bases.c, 10);
        assert_eq!(summary.bases.g, 10);
        assert_eq!(summary.bases.t, 10);
        assert!(n_mask.iter().all(|&byte| byte == 0));
        assert_eq!(summary.qualities.sum_phred, 40 * seq.len() as u64);
    }

    #[test]
    fn fused_pack_falls_back_for_different_quality_len() {
        let mut bases = Vec::new();
        let mut n_mask = Vec::new();
        let summary =
            pack_bases_and_summarize_qualities_into(b"ACGT", b"!I", &mut bases, &mut n_mask)
                .unwrap();

        assert_eq!(bases, vec![0b1110_0100]);
        assert_eq!(n_mask, vec![0]);
        assert_eq!(summary.bases.canonical_bases(), 4);
        assert_eq!(
            summary.qualities,
            QualitySummary {
                len: 2,
                min_phred: Some(0),
                max_phred: Some(40),
                sum_phred: 40,
                q20_bases: 1,
                q30_bases: 1,
            }
        );
    }

    #[test]
    fn rejects_non_printable_quality() {
        let err = summarize_qualities(b"I\nI").unwrap_err();
        assert_eq!(
            err,
            PackError::InvalidQuality {
                offset: 1,
                byte: b'\n',
            }
        );
    }

    #[test]
    fn bins_qualities_by_phred_thresholds() {
        let mut bins = Vec::with_capacity(16);
        let summary = bin_qualities_into(b"!+5?I", &[10, 20, 30], &mut bins).unwrap();
        assert_eq!(bins, vec![0, 1, 2, 3, 3]);
        assert_eq!(summary.len, 5);
        assert_eq!(summary.q20_bases, 3);
        assert!(bins.capacity() >= 16);
    }

    #[test]
    fn binning_slice_checks_output_len() {
        let mut bins = [0; 2];
        let err = bin_qualities_into_slice(b"IIII", &[20, 30], &mut bins).unwrap_err();
        assert_eq!(
            err,
            PackError::OutputTooSmall {
                buffer: PackBuffer::QualityBins,
                needed: 4,
                provided: 2,
            }
        );
    }

    #[test]
    fn rejects_unsorted_thresholds() {
        let mut bins = Vec::new();
        let err = bin_qualities_into(b"IIII", &[20, 10], &mut bins).unwrap_err();
        assert_eq!(err, PackError::UnsortedQualityThresholds { index: 1 });
        assert!(bins.is_empty());
    }

    #[test]
    fn trusted_fastq_exposes_packed_buffers_to_sink() {
        let mut seen = Vec::new();
        pack_trusted_fastq(b"@r0\nACGTN\n+\nIIIII\n", |record| {
            seen.push((
                record.bases.to_vec(),
                record.n_mask.to_vec(),
                record.summary,
            ));
            Ok(())
        })
        .unwrap();

        assert_eq!(seen.len(), 1);
        assert_eq!(seen[0].0, vec![0b1110_0100, 0]);
        assert_eq!(seen[0].1, vec![0b0001_0000]);
        assert_eq!(seen[0].2.bases.n, 1);
    }

    #[test]
    fn trusted_direct_fastq_matches_offset_scan() {
        let input = b"@r0\r\nACGTN\r\n+\r\nIIIII\r\n@r1\nTGCA\n+\n!!!!\n";
        let mut offset = Vec::new();
        let mut direct = Vec::new();

        pack_trusted_fastq(input, |record| {
            offset.push((
                record.name.to_vec(),
                record.bases.to_vec(),
                record.n_mask.to_vec(),
                record.summary,
            ));
            Ok(())
        })
        .unwrap();

        pack_trusted_fastq_direct(input, |record| {
            direct.push((
                record.name.to_vec(),
                record.bases.to_vec(),
                record.n_mask.to_vec(),
                record.summary,
            ));
            Ok(())
        })
        .unwrap();

        assert_eq!(direct, offset);
    }

    #[test]
    fn trusted_direct_stream_handles_slab_carry() {
        let input = b"@r0\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@r1\nNNNN\n+\n!!!!";
        let mut records = 0;
        pack_trusted_fastq_read_direct(
            &input[..],
            FastqConfig {
                slab_size: 16,
                ..FastqConfig::default()
            },
            |record| {
                records += 1;
                assert_eq!(record.summary.bases.len, record.seq.len());
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(records, 2);
    }

    #[test]
    fn trusted_paired_fastq_validates_fast_slash_ids() {
        let r1 = b"@frag/1\nACGT\n+\nIIII\n";
        let r2 = b"@frag/2\nTGCA\n+\nIIII\n";
        let mut pairs = 0;
        pack_trusted_paired_fastq_read(
            &r1[..],
            &r2[..],
            FastqConfig::default(),
            crate::PairValidation::FastSlash,
            |pair| {
                assert_eq!(pair.first.summary.bases.len, 4);
                assert_eq!(pair.second.summary.bases.len, 4);
                pairs += 1;
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(pairs, 1);
    }

    #[test]
    fn trusted_paired_fastq_streams_across_slab_carry() {
        let r1 = b"@frag0/1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@frag1/1\nNNNN\n+\n!!!!\n";
        let r2 = b"@frag0/2\nTGCATGCATGCA\n+\nJJJJJJJJJJJJ\n@frag1/2\nACGT\n+\n####\n";
        let mut pairs = 0;
        pack_trusted_paired_fastq_read(
            &r1[..],
            &r2[..],
            FastqConfig {
                slab_size: 16,
                ..FastqConfig::default()
            },
            crate::PairValidation::FastSlash,
            |pair| {
                assert_eq!(pair.first.summary.bases.len, pair.first.seq.len());
                assert_eq!(pair.second.summary.bases.len, pair.second.seq.len());
                pairs += 1;
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(pairs, 2);
    }

    #[test]
    fn trusted_paired_fastq_rejects_mismatched_ids() {
        let r1 = b"@frag-a/1\nACGT\n+\nIIII\n";
        let r2 = b"@frag-b/2\nTGCA\n+\nIIII\n";
        let err = pack_trusted_paired_fastq_read(
            &r1[..],
            &r2[..],
            FastqConfig::default(),
            crate::PairValidation::FastSlash,
            |_pair| Ok(()),
        )
        .unwrap_err();
        assert!(err.to_string().contains("identifiers do not match"));
    }

    #[test]
    fn trusted_paired_fastq_rejects_different_record_counts() {
        let r1 = b"@frag0/1\nACGT\n+\nIIII\n@frag1/1\nTGCA\n+\nIIII\n";
        let r2 = b"@frag0/2\nTGCA\n+\nIIII\n";
        let err = pack_trusted_paired_fastq_read(
            &r1[..],
            &r2[..],
            FastqConfig {
                slab_size: 16,
                ..FastqConfig::default()
            },
            crate::PairValidation::FastSlash,
            |_pair| Ok(()),
        )
        .unwrap_err();
        assert!(err.to_string().contains("different record counts"));
    }

    #[test]
    fn reports_selected_pack_kernel() {
        #[cfg(feature = "simd")]
        assert!(matches!(
            selected_pack_kernel(),
            PackKernel::PortableSimd | PackKernel::Avx2
        ));
        #[cfg(not(feature = "simd"))]
        assert_eq!(selected_pack_kernel(), PackKernel::Scalar);
    }
}
