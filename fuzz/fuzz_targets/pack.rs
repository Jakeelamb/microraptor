#![no_main]

use libfuzzer_sys::fuzz_target;
use microraptor::pack;

const MAX_INPUT: usize = 256 * 1024;

fuzz_target!(|data: &[u8]| {
    let input = &data[..data.len().min(MAX_INPUT)];
    let (threshold_bytes, payload) = split_thresholds(input);
    let thresholds = sorted_thresholds(threshold_bytes);

    let packed = pack::pack_bases(payload);
    for index in 0..payload.len() {
        let _ = pack::packed_base_at(&packed.bases, &packed.n_mask, index);
    }

    let mut bases = Vec::new();
    let mut n_mask = Vec::new();
    let _ = pack::pack_bases_into(payload, &mut bases, &mut n_mask);

    let mut quality_bins = Vec::new();
    let _ = pack::summarize_qualities(payload);
    let _ = pack::bin_qualities_into(payload, &thresholds, &mut quality_bins);
});

fn split_thresholds(data: &[u8]) -> (&[u8], &[u8]) {
    if data.is_empty() {
        return (&[], &[]);
    }

    let threshold_len = data.first().copied().unwrap_or(0) as usize % 8;
    let split = 1 + threshold_len.min(data.len().saturating_sub(1));
    let (config, payload) = data.split_at(split);
    (&config[1..], payload)
}

fn sorted_thresholds(bytes: &[u8]) -> Vec<u8> {
    let mut thresholds = bytes
        .iter()
        .copied()
        .map(|byte| byte % 94)
        .collect::<Vec<_>>();
    thresholds.sort_unstable();
    thresholds.dedup();
    thresholds
}
