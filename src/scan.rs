#[cfg(all(feature = "simd", target_arch = "x86_64"))]
use std::arch::x86_64::{
    __m256i, _mm256_cmpeq_epi8, _mm256_loadu_si256, _mm256_movemask_epi8, _mm256_set1_epi8,
};

pub(crate) fn scan_newlines(bytes: &[u8], out: &mut Vec<usize>) {
    out.clear();
    scan_newlines_impl(bytes, out);
}

#[cfg(all(feature = "simd", target_arch = "x86_64"))]
fn scan_newlines_impl(bytes: &[u8], out: &mut Vec<usize>) {
    const LANES: usize = 32;
    let mut i = 0;
    if std::is_x86_feature_detected!("avx2") {
        unsafe {
            let needle = _mm256_set1_epi8(b'\n' as i8);
            while i + LANES <= bytes.len() {
                let chunk = _mm256_loadu_si256(bytes.as_ptr().add(i).cast::<__m256i>());
                let mut mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(chunk, needle)) as u32;
                while mask != 0 {
                    let bit = mask.trailing_zeros() as usize;
                    out.push(i + bit);
                    mask &= mask - 1;
                }
                i += LANES;
            }
        }
    }
    out.extend(memchr::memchr_iter(b'\n', &bytes[i..]).map(|offset| i + offset));
}

#[cfg(not(all(feature = "simd", target_arch = "x86_64")))]
fn scan_newlines_impl(bytes: &[u8], out: &mut Vec<usize>) {
    out.extend(memchr::memchr_iter(b'\n', bytes));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scans_newlines() {
        let mut out = Vec::new();
        scan_newlines(b"a\nbc\n\nz", &mut out);
        assert_eq!(out, vec![1, 4, 5]);
    }
}
