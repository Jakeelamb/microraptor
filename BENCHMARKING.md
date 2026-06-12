# Benchmarking

Microraptor has three benchmark surfaces:

- `cargo bench --all-features`: nightly microbenchmarks using Rust's built-in
  benchmark harness.
- `cargo run --release --bin microraptor-bench -- ...`: release-mode throughput
  benchmark with table or JSON output.
- `scripts/profile-perf.sh`: Linux `perf stat` plus sampled call graph output.

The benchmark binary generates deterministic synthetic FASTQ in memory, then
measures the same parser and side-channel APIs used by downstream crates. It
reports best-of-N wall time to reduce noise from scheduler spikes.

For real input files, pass `--input PATH` or set `MICRORAPTOR_INPUT`. The file
path goes through `open_fastq_with_config`, so raw FASTQ, gzip FASTQ, and BGZF
FASTQ use the same auto-detection path as library callers.

## Fast Commands

```bash
cargo test --all
cargo clippy --all-targets --all-features -- -D warnings
cargo bench --all-features
scripts/bench.sh
```

For machine-readable output:

```bash
cargo run --release --bin microraptor-bench -- --records 500000 --iters 7 --json
```

For a real dataset:

```bash
cargo run --release --bin microraptor-bench -- --input reads.fastq.gz --iters 5
MICRORAPTOR_INPUT=reads.fastq.gz scripts/bench.sh
```

## Profiling

```bash
scripts/profile-perf.sh
```

Outputs:

- `target/profiles/microraptor-bench.perf.data`
- `target/profiles/microraptor-bench.perf.txt`

Use larger inputs when looking for stable instruction-cache, branch, and memory
behavior:

```bash
MICRORAPTOR_RECORDS=2000000 MICRORAPTOR_ITERS=3 scripts/profile-perf.sh
```

Profile a real input file with the same perf commands:

```bash
MICRORAPTOR_INPUT=reads.fastq.gz MICRORAPTOR_ITERS=3 scripts/profile-perf.sh
```

## Interpreting Results

The benchmark table reports:

- `input_mib_s`: throughput over the input representation being measured. For
  compressed inputs this is compressed bytes per second.
- `records_s`: FASTQ records consumed per second.
- `bases_s`: sequence bases consumed per second.
- `checksum`: optimization guard; changes indicate behavior changed.

Do not compare raw and compressed `input_mib_s` directly. For biological
pipeline planning, `records_s` and `bases_s` are the more useful common units.

The `bgzf-parallel` row uses the bounded streaming `BgzfParallelReader`, not the
older whole-input decompression helper.
