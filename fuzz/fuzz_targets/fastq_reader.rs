#![no_main]

use libfuzzer_sys::fuzz_target;
use microraptor::{FastqConfig, FastqReader, pack};

const MAX_BATCHES: usize = 256;

fuzz_target!(|data: &[u8]| {
    let (config_bytes, fastq_bytes) = split_config(data);
    let mut reader = FastqReader::with_config(
        fastq_bytes,
        FastqConfig {
            slab_size: slab_size(config_bytes),
            validate: validate(config_bytes),
            ..FastqConfig::default()
        },
    );

    for _ in 0..MAX_BATCHES {
        let Ok(Some(batch)) = reader.next_batch() else {
            return;
        };

        for record in batch.records() {
            let seq = record.seq();
            let qual = record.qual();

            let _ = record.name();
            let _ = record.name_without_at();
            let _ = record.id_token();
            let _ = record.pair_normalized_id();
            let _ = record.plus();
            let _ = pack::pack_bases(seq);
            let _ = pack::summarize_qualities(qual);
        }
    }
});

fn split_config(data: &[u8]) -> (&[u8], &[u8]) {
    if data.len() < 4 {
        (data, &[])
    } else {
        data.split_at(4)
    }
}

fn slab_size(config: &[u8]) -> usize {
    let low = config.first().copied().unwrap_or(0) as usize;
    let high = config.get(1).copied().unwrap_or(0) as usize;
    1024 + ((high << 8) | low) % (64 * 1024)
}

fn validate(config: &[u8]) -> bool {
    config.get(2).copied().unwrap_or(1) & 1 == 1
}
