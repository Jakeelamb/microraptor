use std::fmt;
use std::io::Read;
#[cfg(feature = "simd")]
use std::simd::{Simd, cmp::SimdPartialOrd};

use crate::scan::scan_newlines;
use crate::{FastqConfig, FastqError, FastqPosition, Result as FastqResult};

/// Summary of a packed DNA sequence.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BaseSummary {
    pub len: usize,
    pub a: usize,
    pub c: usize,
    pub g: usize,
    pub t: usize,
    pub n: usize,
}

impl BaseSummary {
    pub fn gc_bases(self) -> usize {
        self.c + self.g
    }

    pub fn canonical_bases(self) -> usize {
        self.a + self.c + self.g + self.t
    }

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
    pub bases: Vec<u8>,
    pub n_mask: Vec<u8>,
    pub summary: BaseSummary,
}

impl PackedSequence {
    pub fn len(&self) -> usize {
        self.summary.len
    }

    pub fn is_empty(&self) -> bool {
        self.summary.is_empty()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackedBase {
    A,
    C,
    G,
    T,
    N,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackBuffer {
    Bases,
    NMask,
    QualityBins,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackError {
    OutputTooSmall {
        buffer: PackBuffer,
        needed: usize,
        provided: usize,
    },
    InvalidQuality {
        offset: usize,
        byte: u8,
    },
    UnsortedQualityThresholds {
        index: usize,
    },
    TooManyQualityThresholds {
        count: usize,
    },
}

/// Phred+33 quality summary.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct QualitySummary {
    pub len: usize,
    pub min_phred: Option<u8>,
    pub max_phred: Option<u8>,
    pub sum_phred: u64,
    pub q20_bases: usize,
    pub q30_bases: usize,
}

impl QualitySummary {
    pub fn mean_phred(self) -> Option<f64> {
        if self.len == 0 {
            None
        } else {
            Some(self.sum_phred as f64 / self.len as f64)
        }
    }

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

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PackedRecordSummary {
    pub bases: BaseSummary,
    pub qualities: QualitySummary,
}

#[derive(Debug, Clone, Copy)]
pub struct TrustedPackedRecord<'a> {
    pub name: &'a [u8],
    pub seq: &'a [u8],
    pub qual: &'a [u8],
    pub bases: &'a [u8],
    pub n_mask: &'a [u8],
    pub summary: PackedRecordSummary,
}

#[derive(Debug, Clone, Copy)]
pub struct TrustedPackedPair<'a> {
    pub first: TrustedPackedRecord<'a>,
    pub second: TrustedPackedRecord<'a>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackKernel {
    Scalar,
    PortableSimd,
}

pub fn selected_pack_kernel() -> PackKernel {
    select_pack_kernel()
}

#[cfg(feature = "simd")]
fn select_pack_kernel() -> PackKernel {
    PackKernel::PortableSimd
}

#[cfg(not(feature = "simd"))]
fn select_pack_kernel() -> PackKernel {
    PackKernel::Scalar
}

#[derive(Debug, Clone, Copy)]
pub struct TrustedPackSlab {
    pub records: u64,
}

pub trait TrustedPackSink {
    fn record(&mut self, record: TrustedPackedRecord<'_>) -> FastqResult<()>;

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

pub fn pack_trusted_fastq(
    input: &[u8],
    on_record: impl FnMut(TrustedPackedRecord<'_>) -> FastqResult<()>,
) -> FastqResult<()> {
    pack_trusted_fastq_sink(input, on_record)
}

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

pub fn pack_trusted_fastq_read<R: Read>(
    mut reader: R,
    config: FastqConfig,
    on_record: impl FnMut(TrustedPackedRecord<'_>) -> FastqResult<()>,
) -> FastqResult<()> {
    pack_trusted_fastq_read_sink(&mut reader, config, on_record)
}

pub fn pack_trusted_fastq_read_sink<R: Read>(
    mut reader: R,
    config: FastqConfig,
    mut sink: impl TrustedPackSink,
) -> FastqResult<()> {
    let slab_size = config.slab_size.max(1024);
    let mut buf = vec![0_u8; slab_size];
    let mut len = 0;
    let mut eof = false;
    let mut base_offset = 0_u64;
    let mut record_index = 0_u64;
    let mut newlines = Vec::with_capacity(slab_size / 48);
    let mut bases = Vec::new();
    let mut n_mask = Vec::new();

    loop {
        while !eof && len < slab_size {
            let n = reader.read(&mut buf[len..slab_size])?;
            if n == 0 {
                eof = true;
                break;
            }
            len += n;
        }

        let slab = pack_trusted_fastq_slab(
            &buf[..len],
            SlabContext {
                base_offset,
                first_record_index: record_index,
                eof,
            },
            &mut newlines,
            &mut bases,
            &mut n_mask,
            &mut sink,
        )?;
        if slab.records != 0 {
            sink.slab(TrustedPackSlab {
                records: slab.records,
            })?;
        }
        record_index += slab.records;

        if slab.next_start == len {
            base_offset += len as u64;
            len = 0;
        } else {
            let carry = len - slab.next_start;
            if slab.next_start == 0 && carry == slab_size && !eof {
                return Err(FastqError::RecordTooLarge { slab_size });
            }
            buf.copy_within(slab.next_start..len, 0);
            base_offset += slab.next_start as u64;
            len = carry;
        }

        if eof {
            if len == 0 {
                return Ok(());
            }
            return Err(FastqError::RecordTooLarge { slab_size });
        }
    }
}

pub fn pack_trusted_paired_fastq_read<R1: Read, R2: Read>(
    first: R1,
    second: R2,
    config: FastqConfig,
    pair_validation: crate::PairValidation,
    mut on_pair: impl FnMut(TrustedPackedPair<'_>) -> FastqResult<()>,
) -> FastqResult<()> {
    let mut first_records = Vec::new();
    let mut second_records = Vec::new();
    pack_trusted_fastq_read(first, config.clone(), |record| {
        first_records.push(OwnedPackedRecord::from(record));
        Ok(())
    })?;
    pack_trusted_fastq_read(second, config, |record| {
        second_records.push(OwnedPackedRecord::from(record));
        Ok(())
    })?;

    if first_records.len() != second_records.len() {
        return Err(FastqError::Format(
            "paired FASTQ inputs have different record counts".into(),
        ));
    }

    for (index, (first, second)) in first_records.iter().zip(&second_records).enumerate() {
        if pair_validation != crate::PairValidation::None
            && !trusted_pair_ids_match(&first.name, &second.name, pair_validation)
        {
            return Err(FastqError::FormatAt {
                message: "paired FASTQ record identifiers do not match".into(),
                position: FastqPosition::new(0, index as u64, 0),
            });
        }
        on_pair(TrustedPackedPair {
            first: first.as_borrowed(),
            second: second.as_borrowed(),
        })?;
    }

    Ok(())
}

pub const fn packed_base_len(base_count: usize) -> usize {
    base_count / 4 + if base_count.is_multiple_of(4) { 0 } else { 1 }
}

pub const fn bit_mask_len(bit_count: usize) -> usize {
    bit_count / 8 + if bit_count.is_multiple_of(8) { 0 } else { 1 }
}

const BASE_N: u8 = 4;
const BASE_LUT: [u8; 256] = base_lut();

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

pub fn pack_bases_into(seq: &[u8], bases: &mut Vec<u8>, n_mask: &mut Vec<u8>) -> BaseSummary {
    bases.clear();
    n_mask.clear();
    bases.resize(packed_base_len(seq.len()), 0);
    n_mask.resize(bit_mask_len(seq.len()), 0);
    pack_bases_exact(seq, bases, n_mask)
}

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

pub fn is_masked(n_mask: &[u8], index: usize) -> Option<bool> {
    let mask_byte = *n_mask.get(index / 8)?;
    Some(((mask_byte >> (index % 8)) & 1) != 0)
}

pub fn summarize_qualities(qualities: &[u8]) -> Result<QualitySummary, PackError> {
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

#[derive(Debug, Clone)]
struct OwnedPackedRecord {
    name: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    bases: Vec<u8>,
    n_mask: Vec<u8>,
    summary: PackedRecordSummary,
}

impl OwnedPackedRecord {
    fn from(record: TrustedPackedRecord<'_>) -> Self {
        Self {
            name: record.name.to_vec(),
            seq: record.seq.to_vec(),
            qual: record.qual.to_vec(),
            bases: record.bases.to_vec(),
            n_mask: record.n_mask.to_vec(),
            summary: record.summary,
        }
    }

    fn as_borrowed(&self) -> TrustedPackedRecord<'_> {
        TrustedPackedRecord {
            name: &self.name,
            seq: &self.seq,
            qual: &self.qual,
            bases: &self.bases,
            n_mask: &self.n_mask,
            summary: self.summary,
        }
    }
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
    while chunk_index < full_chunks {
        let c0 = BASE_LUT[usize::from(seq[base_index])];
        let c1 = BASE_LUT[usize::from(seq[base_index + 1])];
        let c2 = BASE_LUT[usize::from(seq[base_index + 2])];
        let c3 = BASE_LUT[usize::from(seq[base_index + 3])];

        let mut packed = 0_u8;
        packed |= pack_code(c0, base_index, 0, &mut summary, n_mask);
        packed |= pack_code(c1, base_index, 1, &mut summary, n_mask);
        packed |= pack_code(c2, base_index, 2, &mut summary, n_mask);
        packed |= pack_code(c3, base_index, 3, &mut summary, n_mask);
        bases[chunk_index] = packed;
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

    let mut bases_summary = BaseSummary {
        len: seq.len(),
        ..BaseSummary::default()
    };

    let full_chunks = seq.len() / 4;
    let mut chunk_index = 0;
    let mut base_index = 0;
    while chunk_index < full_chunks {
        let c0 = BASE_LUT[usize::from(seq[base_index])];
        let c1 = BASE_LUT[usize::from(seq[base_index + 1])];
        let c2 = BASE_LUT[usize::from(seq[base_index + 2])];
        let c3 = BASE_LUT[usize::from(seq[base_index + 3])];

        let mut packed = 0_u8;
        packed |= pack_code(c0, base_index, 0, &mut bases_summary, n_mask);
        packed |= pack_code(c1, base_index, 1, &mut bases_summary, n_mask);
        packed |= pack_code(c2, base_index, 2, &mut bases_summary, n_mask);
        packed |= pack_code(c3, base_index, 3, &mut bases_summary, n_mask);
        bases[chunk_index] = packed;
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
        index += 1;
    }

    Ok(PackedRecordSummary {
        bases: bases_summary,
        qualities: summarize_qualities(qualities)?,
    })
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
fn pack_code(
    code: u8,
    base_index: usize,
    offset: usize,
    summary: &mut BaseSummary,
    n_mask: &mut [u8],
) -> u8 {
    if code < BASE_N {
        add_base_count(summary, code);
        code << (offset * 2)
    } else {
        summary.n += 1;
        let index = base_index + offset;
        n_mask[index / 8] |= 1 << (index % 8);
        0
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
    fn reports_selected_pack_kernel() {
        #[cfg(feature = "simd")]
        assert_eq!(selected_pack_kernel(), PackKernel::PortableSimd);
        #[cfg(not(feature = "simd"))]
        assert_eq!(selected_pack_kernel(), PackKernel::Scalar);
    }
}
