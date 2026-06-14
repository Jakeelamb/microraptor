#![no_main]

use libfuzzer_sys::fuzz_target;
use microraptor::{
    FastaConfig, FastaReader, build_fasta_index, count_fasta_bytes, count_fasta_read,
    detect_fasta_shape, visit_fasta_bytes, visit_fasta_bytes_auto,
};

const MAX_BATCHES: usize = 256;
const MAX_INPUT: usize = 256 * 1024;

fuzz_target!(|data: &[u8]| {
    let input = &data[..data.len().min(MAX_INPUT)];
    let (config_bytes, fasta_bytes) = split_config(input);
    let mut reader = FastaReader::with_config(
        fasta_bytes,
        FastaConfig {
            batch_records: batch_records(config_bytes),
            buffer_size: buffer_size(config_bytes),
            expected_seq_len: expected_seq_len(config_bytes),
        },
    );

    for _ in 0..MAX_BATCHES {
        let Ok(Some(batch)) = reader.next_batch() else {
            break;
        };
        for record in batch.records() {
            let _ = record.name();
            let _ = record.name_without_gt();
            let _ = record.id_token();
            let _ = record.seq();
        }
    }

    let _ = count_fasta_read(fasta_bytes);
    let _ = count_fasta_bytes(fasta_bytes);
    let _ = detect_fasta_shape(fasta_bytes);
    let _ = visit_fasta_bytes(fasta_bytes, |_| Ok(()));
    let _ = visit_fasta_bytes_auto(fasta_bytes, |_| Ok(()));
    let _ = build_fasta_index(fasta_bytes);
});

fn split_config(data: &[u8]) -> (&[u8], &[u8]) {
    if data.len() < 4 {
        (data, &[])
    } else {
        data.split_at(4)
    }
}

fn batch_records(config: &[u8]) -> usize {
    1 + config.first().copied().unwrap_or(0) as usize % 128
}

fn buffer_size(config: &[u8]) -> usize {
    let low = config.get(1).copied().unwrap_or(0) as usize;
    let high = config.get(2).copied().unwrap_or(0) as usize;
    1024 + ((high << 8) | low) % (64 * 1024)
}

fn expected_seq_len(config: &[u8]) -> usize {
    config.get(3).copied().unwrap_or(0) as usize * 16
}
