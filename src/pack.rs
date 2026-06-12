use std::fmt;

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
    let mut summary = QualityAccumulator::default();
    for (offset, &byte) in qualities.iter().enumerate() {
        summary.observe(byte, offset)?;
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
    let mut quality = QualityAccumulator::default();

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
        quality.observe(qualities[base_index], base_index)?;
        quality.observe(qualities[base_index + 1], base_index + 1)?;
        quality.observe(qualities[base_index + 2], base_index + 2)?;
        quality.observe(qualities[base_index + 3], base_index + 3)?;
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
        quality.observe(qualities[index], index)?;
        index += 1;
    }

    Ok(PackedRecordSummary {
        bases: bases_summary,
        qualities: quality.finish(),
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
}
