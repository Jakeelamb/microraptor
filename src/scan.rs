#[cfg(feature = "simd")]
use std::simd::{Simd, cmp::SimdPartialEq};

pub(crate) fn scan_newlines(bytes: &[u8], out: &mut Vec<usize>) {
    out.clear();
    scan_newlines_impl(bytes, out);
}

#[cfg(feature = "simd")]
fn scan_newlines_impl(bytes: &[u8], out: &mut Vec<usize>) {
    const LANES: usize = 64;
    type Chunk = Simd<u8, LANES>;

    let needle = Chunk::splat(b'\n');
    let mut i = 0;
    while i + LANES <= bytes.len() {
        let chunk = Chunk::from_slice(&bytes[i..i + LANES]);
        let mut mask = chunk.simd_eq(needle).to_bitmask();
        while mask != 0 {
            let bit = mask.trailing_zeros() as usize;
            out.push(i + bit);
            mask &= mask - 1;
        }
        i += LANES;
    }
    scan_newlines_scalar_offset(&bytes[i..], i, out);
}

#[cfg(not(feature = "simd"))]
fn scan_newlines_impl(bytes: &[u8], out: &mut Vec<usize>) {
    out.extend(memchr::memchr_iter(b'\n', bytes));
}

#[cfg(feature = "simd")]
fn scan_newlines_scalar_offset(bytes: &[u8], offset: usize, out: &mut Vec<usize>) {
    for (i, &b) in bytes.iter().enumerate() {
        if b == b'\n' {
            out.push(offset + i);
        }
    }
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
