#![feature(test)]

extern crate test;

use microraptor::benchutil::{consume_fastq, synthetic_fastq};
use microraptor::pack::pack_bases_and_summarize_qualities_into;
use microraptor::{FastqConfig, FastqReader};
use test::{Bencher, black_box};

const RECORDS: usize = 50_000;
const READ_LEN: usize = 150;

#[bench]
fn parse_raw_strict_8m_slab(b: &mut Bencher) {
    let input = synthetic_fastq(RECORDS, READ_LEN);
    b.bytes = input.len() as u64;
    b.iter(|| {
        let mut reader = FastqReader::with_config(
            &input[..],
            FastqConfig {
                slab_size: 8 * 1024 * 1024,
                validate: true,
                ..FastqConfig::default()
            },
        );
        black_box(consume_fastq(&mut reader).unwrap());
    });
}

#[bench]
fn parse_raw_no_validate_8m_slab(b: &mut Bencher) {
    let input = synthetic_fastq(RECORDS, READ_LEN);
    b.bytes = input.len() as u64;
    b.iter(|| {
        let mut reader = FastqReader::with_config(
            &input[..],
            FastqConfig {
                slab_size: 8 * 1024 * 1024,
                validate: false,
                ..FastqConfig::default()
            },
        );
        black_box(consume_fastq(&mut reader).unwrap());
    });
}

#[bench]
fn parse_and_pack_seq_qual(b: &mut Bencher) {
    let input = synthetic_fastq(RECORDS, READ_LEN);
    b.bytes = input.len() as u64;
    b.iter(|| {
        let mut reader = FastqReader::new(&input[..]);
        let mut packed = Vec::new();
        let mut mask = Vec::new();
        let mut checksum = 0_u64;
        while let Some(batch) = reader.next_batch().unwrap() {
            for record in batch.records() {
                let summary = pack_bases_and_summarize_qualities_into(
                    record.seq(),
                    record.qual(),
                    &mut packed,
                    &mut mask,
                )
                .unwrap();
                checksum = checksum
                    .wrapping_add(summary.bases.canonical_bases() as u64)
                    .wrapping_add(summary.qualities.sum_phred);
            }
        }
        black_box(checksum);
    });
}

#[cfg(feature = "bgzf")]
#[bench]
fn bgzf_parallel_decode_then_parse(b: &mut Bencher) {
    let input = synthetic_fastq(RECORDS, READ_LEN);
    let encoded = microraptor::compress_bgzf_parallel(&input, 4).unwrap();
    b.bytes = input.len() as u64;
    b.iter(|| {
        let decoded = microraptor::decompress_bgzf_parallel(&encoded[..], 4).unwrap();
        let mut reader = FastqReader::new(std::io::Cursor::new(decoded));
        black_box(consume_fastq(&mut reader).unwrap());
    });
}
