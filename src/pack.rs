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

use crate::fastq_frame::{self, RecordLines, RecordValidation};
use crate::scan::scan_newlines;
use crate::{FastqConfig, FastqError, Result as FastqResult};

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
    /// Header line, including the leading `@` and excluding the newline.
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
    /// Legacy portable-SIMD implementation marker.
    ///
    /// Current `simd` builds use stable `std::arch` paths where available and
    /// otherwise fall back to [`Scalar`](Self::Scalar).
    PortableSimd,
    /// x86-64 AVX2 implementation.
    Avx2,
}

/// Return the base-packing kernel selected for this build and host.
pub fn selected_pack_kernel() -> PackKernel {
    select_pack_kernel()
}

#[cfg(all(feature = "simd", target_arch = "x86_64"))]
fn select_pack_kernel() -> PackKernel {
    if std::is_x86_feature_detected!("avx2") {
        return PackKernel::Avx2;
    }
    PackKernel::Scalar
}

#[cfg(not(all(feature = "simd", target_arch = "x86_64")))]
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

            fastq_frame::validate_record(
                first,
                first_reader.base_offset,
                first_record_index,
                RecordValidation::TRUSTED_PACK,
            )?;
            fastq_frame::validate_record(
                second,
                second_reader.base_offset,
                second_record_index,
                RecordValidation::TRUSTED_PACK,
            )?;

            if pair_validation != crate::PairValidation::None
                && !trusted_pair_ids_match(first.name.bytes, second.name.bytes, pair_validation)
            {
                return Err(fastq_frame::format_at(
                    "paired FASTQ record identifiers do not match",
                    0,
                    0,
                    pair_index,
                    0,
                ));
            }

            let first_summary = pack_bases_and_summarize_qualities_into(
                first.seq.bytes,
                first.qual.bytes,
                &mut first_bases,
                &mut first_n_mask,
            )
            .map_err(|err| {
                fastq_frame::format_at(
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
                fastq_frame::format_at(
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
            return Err(fastq_frame::format_at(
                "truncated FASTQ record",
                self.base_offset,
                fastq_frame::line_start(&self.newlines, (self.line_count / 4) * 4),
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

    fn record_lines(&self, index: usize) -> RecordLines<'_> {
        let line = index * 4;
        fastq_frame::record_lines(&self.buf[..self.len], &self.newlines, line)
    }

    fn consume_records(&mut self, records: usize) -> FastqResult<()> {
        if records == 0 {
            return Ok(());
        }
        let consumed_lines = records * 4;
        let next_start = if self.eof && consumed_lines == self.line_count {
            self.len
        } else {
            fastq_frame::line_start(&self.newlines, consumed_lines)
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
    let summary = pack_bases_exact_zeroed(seq, &mut bases, &mut n_mask);
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
    pack_bases_exact_zeroed(seq, bases, n_mask)
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
        pack_bases_and_qualities_exact_zeroed(seq, qualities, bases, n_mask)
    } else {
        Ok(PackedRecordSummary {
            bases: pack_bases_exact_zeroed(seq, bases, n_mask),
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

    let layout = fastq_frame::slab_line_layout(input, newlines, context.eof);
    if context.eof && !layout.line_count.is_multiple_of(4) {
        let record_index = context.first_record_index + (layout.line_count / 4) as u64;
        return Err(fastq_frame::format_at(
            "truncated FASTQ record",
            context.base_offset,
            fastq_frame::line_start(newlines, (layout.line_count / 4) * 4),
            record_index,
            (layout.line_count % 4) as u8,
        ));
    }

    for line in (0..layout.complete_lines).step_by(4) {
        let record_index = context.first_record_index + (line / 4) as u64;
        let record = fastq_frame::record_lines(input, newlines, line);

        observe_trusted_packed_record(
            record,
            context.base_offset,
            record_index,
            bases,
            n_mask,
            sink,
        )?;
        records += 1;
    }

    if layout.complete_lines == layout.line_count && context.eof {
        Ok(SlabResult {
            next_start: input.len(),
            records,
        })
    } else {
        let next_start = fastq_frame::line_start(newlines, layout.complete_lines);
        if layout.complete_lines == layout.line_count && next_start == input.len() {
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
        let Some(name) = fastq_frame::direct_line(input, &mut cursor, context.eof) else {
            return Ok(SlabResult {
                next_start: record_start,
                records,
            });
        };
        let Some(seq) = fastq_frame::direct_line(input, &mut cursor, context.eof) else {
            return incomplete_or_truncated_direct(input, context, record_start, records, 1);
        };
        let Some(plus) = fastq_frame::direct_line(input, &mut cursor, context.eof) else {
            return incomplete_or_truncated_direct(input, context, record_start, records, 2);
        };
        let Some(qual) = fastq_frame::direct_line(input, &mut cursor, context.eof) else {
            return incomplete_or_truncated_direct(input, context, record_start, records, 3);
        };
        let record = RecordLines {
            name,
            seq,
            plus,
            qual,
        };

        observe_trusted_packed_record(
            record,
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
        Err(fastq_frame::format_at(
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

fn observe_trusted_packed_record(
    record: RecordLines<'_>,
    base_offset: u64,
    record_index: u64,
    bases: &mut Vec<u8>,
    n_mask: &mut Vec<u8>,
    sink: &mut impl TrustedPackSink,
) -> FastqResult<()> {
    fastq_frame::validate_record(
        record,
        base_offset,
        record_index,
        RecordValidation::TRUSTED_PACK,
    )?;

    let summary =
        pack_bases_and_summarize_qualities_into(record.seq.bytes, record.qual.bytes, bases, n_mask)
            .map_err(|err| {
                fastq_frame::format_at(
                    err.to_string(),
                    base_offset,
                    record.qual.start,
                    record_index,
                    3,
                )
            })?;
    sink.record(TrustedPackedRecord {
        name: record.name.bytes,
        seq: record.seq.bytes,
        qual: record.qual.bytes,
        bases: &bases[..],
        n_mask: &n_mask[..],
        summary,
    })
}

fn trusted_pair_ids_match(
    first_name: &[u8],
    second_name: &[u8],
    mode: crate::PairValidation,
) -> bool {
    match mode {
        crate::PairValidation::None => true,
        crate::PairValidation::FastSlash => {
            fastq_frame::fast_slash_pair_ids_match(first_name, second_name).unwrap_or_else(|| {
                fastq_frame::normalized_pair_id(first_name)
                    == fastq_frame::normalized_pair_id(second_name)
            })
        }
        crate::PairValidation::Full => {
            fastq_frame::normalized_pair_id(first_name)
                == fastq_frame::normalized_pair_id(second_name)
        }
    }
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

    pack_bases_exact_zeroed(seq, bases, n_mask)
}

fn pack_bases_exact_zeroed(seq: &[u8], bases: &mut [u8], n_mask: &mut [u8]) -> BaseSummary {
    let bases_needed = packed_base_len(seq.len());
    if bases_needed == 0 {
        return BaseSummary::default();
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

fn pack_bases_and_qualities_exact_zeroed(
    seq: &[u8],
    qualities: &[u8],
    bases: &mut [u8],
    n_mask: &mut [u8],
) -> Result<PackedRecordSummary, PackError> {
    let bases_needed = packed_base_len(seq.len());
    if bases_needed == 0 {
        return Ok(PackedRecordSummary::default());
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
    let mut codes = [BASE_N; 16];
    let mut i = 0;
    while i < 16 {
        codes[i] = BASE_LUT[usize::from(seq[i])];
        i += 1;
    }
    codes
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

#[cfg(test)]
mod tests;
