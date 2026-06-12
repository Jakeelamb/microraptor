# Benchmarking

Microraptor has five benchmark surfaces:

- `cargo bench --all-features`: nightly microbenchmarks using Rust's built-in
  benchmark harness.
- `cargo run --release --bin microraptor-bench -- ...`: release-mode throughput
  benchmark with table or JSON output.
- `scripts/profile-perf.sh`: Linux `perf stat` plus sampled call graph output.
- `scripts/benchmark-gauntlet.sh`: generated real-file gauntlet covering raw,
  gzip, BGZF, paired R1/R2, and interleaved inputs.
- `scripts/profile-hotpath.sh`: isolated parse-vs-pack profiling output.

The benchmark binary generates deterministic synthetic FASTQ in memory, then
measures the same parser and side-channel APIs used by downstream crates. It
reports best-of-N wall time to reduce noise from scheduler spikes.

For real input files, pass `--input PATH` or set `MICRORAPTOR_INPUT`. The file
path goes through `open_fastq_with_config`, so raw FASTQ, gzip FASTQ, and BGZF
FASTQ use the same auto-detection path as library callers.

The FASTQ parser expects the common four-line record shape: name, sequence,
plus, quality. Multiline sequence or quality fields are not supported.

Paired input coverage is ordered-pair coverage: generated R1/R2 files with
matching order and interleaved files with adjacent mates. The stateful paired
reader handles different batch boundaries and validates normalized IDs, but it
does not synchronize reordered mates.

## Fast Commands

```bash
cargo test --all
cargo clippy --all-targets --all-features -- -D warnings
cargo clippy --all-targets --no-default-features -- -D warnings
cargo bench --all-features
scripts/bench.sh
scripts/benchmark-gauntlet.sh
```

For machine-readable output:

```bash
cargo run --release --bin microraptor-bench -- --records 500000 --iters 7 --json
cargo run --release --bin microraptor-bench -- --records 500000 --mode parse --json
cargo run --release --bin microraptor-bench -- --records 500000 --mode pack --json
```

Synthetic `--mode pack` reports both the safe parser-backed `pack-seq-qual`
row and the `trusted-pack-seq-qual` row. The trusted row is a narrow raw
four-line FASTQ path that reuses the SIMD newline scanner but skips batch
record construction before packing.

For a real dataset:

```bash
cargo run --release --bin microraptor-bench -- --input reads.fastq.gz --iters 5
cargo run --release --bin microraptor-bench -- --paired-inputs r1.fastq.gz r2.fastq.gz --iters 5
MICRORAPTOR_INPUT=reads.fastq.gz scripts/bench.sh
MICRORAPTOR_INPUT=reads.fastq.gz MICRORAPTOR_MODE=parse scripts/bench.sh
```

The gauntlet writes:

- `target/bench-results/microraptor-gauntlet.jsonl`
- `target/bench-results/microraptor-gauntlet.md`

It also records whether optional external comparators such as `seqkit` or
`fastp` were installed and runnable.

Build with `--features libdeflate` or `--all-features` to include explicit
`bgzf-libdeflate-*` rows for synthetic BGZF and `file-bgzf-libdeflate-*` rows
for real `.bgz` inputs. Ordinary gzip remains on the streaming flate2 path;
`open_fastq_gzip_libdeflate` is an explicit buffered path for bounded gzip
inputs. When the `libdeflate` feature is enabled, the normal BGZF auto-open path
uses libdeflate inflate by default; use `open_fastq_bgzf_flate2` or
`open_fastq_bgzf_with_backend` when comparing or forcing a backend. BGZF output
can use libdeflate through `BgzfDeflateBackend`. Use `open_fastq_bgzf_adaptive`
when you want the crate to keep small BGZF inputs serial and switch to the
bounded parallel reader only past the built-in size threshold.

`microraptor-bench --paired-inputs` uses typed file openers and
`PairValidation::FastSlash` for ordered `/1` and `/2` mate IDs. That benchmark
path is meant to represent the high-performance internal pipeline mode, not the
most defensive public opener configuration.

`build_bgzf_index` records compressed block offsets and uncompressed block
starts. Use `virtual_offset_for_uncompressed_offset` to plan a seek, then
`BgzfSeekReader::seek_virtual_offset` to resume reading from that BGZF virtual
offset.

## Profiling

```bash
scripts/profile-perf.sh
scripts/profile-hotpath.sh
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

For isolated parse/pack hot-path evidence:

```bash
MICRORAPTOR_PROFILE_RECORDS=1000000 MICRORAPTOR_PROFILE_ITERS=3 scripts/profile-hotpath.sh
MICRORAPTOR_PROFILE_INPUT=reads.fastq.gz scripts/profile-hotpath.sh
```

`profile-hotpath.sh` writes:

- `target/profiles/microraptor-hotpath.jsonl`
- `target/profiles/microraptor-parse.perf-stat.txt` when `perf` is permitted
- `target/profiles/microraptor-pack.perf-stat.txt` when `perf` is permitted

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

## CI Parity

The GitHub Actions workflow uses nightly Rust from `rust-toolchain.toml` and
checks:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo clippy --all-targets --no-default-features -- -D warnings
cargo test --all-features
cargo test --no-default-features
cargo fuzz build
```

`cargo fuzz build` only compiles the fuzz targets. It does not run long fuzzing
campaigns in CI.
