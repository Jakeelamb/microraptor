# Benchmarking

Microraptor has fifteen benchmark surfaces:

- `cargo +nightly bench --all-features`: nightly microbenchmarks using Rust's built-in
  benchmark harness.
- `cargo run --release --bin microraptor-bench -- ...`: release-mode throughput
  benchmark with table or JSON output.
- `scripts/profile-perf.sh`: Linux `perf stat` plus sampled call graph output.
- `scripts/benchmark-gauntlet.sh`: generated real-file gauntlet covering raw,
  gzip, BGZF, paired R1/R2, and interleaved inputs.
- `scripts/render-benchmark-report.sh`: JSONL-to-Markdown/SVG renderer for
  publishable benchmark snapshots.
- `scripts/check-benchmark-snapshots.sh`: checked-artifact verifier that
  re-renders stored gauntlet JSONL and Rust/FASTA peer TSV snapshots, validates
  FASTA gauntlet artifacts, and diffs summaries/figures.
- `scripts/benchmark-common.sh`: shared helper surface for generated benchmark
  harnesses, including thread caps, sanitized artifact copies, display paths,
  and feature dependency rendering.
- `scripts/benchmark-rust-peers.sh`: script-generated Rust parser-library peer
  comparison against `seq_io`, `noodles-fastq`, and `bio`.
- `scripts/benchmark-fasta-peers.sh`: script-generated FASTA parser-library
  peer comparison against `seq_io` and `bio`, including strict two-line FASTA
  fast paths.
- `scripts/check-replication-host.sh`: local host/toolchain/comparator
  preflight for release and benchmark regeneration.
- `scripts/discover-local-benchmark-corpus.sh`: local biological FASTQ corpus
  discovery for the `~/Projects/Benchmarks` workspace.
- `scripts/profile-hotpath.sh`: isolated parse-vs-pack profiling output.
- `scripts/benchmark-rust-peer-size-sweep.sh`: increasing-input-size Rust parser
  framework sweep for publication-style time-vs-size figures.
- `scripts/benchmark-fasta-peer-size-sweep.sh`: increasing-input-size FASTA
  parser framework sweep for raw and gzip two-line FASTA figures.
- `scripts/benchmark-fasta-gauntlet.sh`: FASTA shape/transport gauntlet covering
  two-line DNA, wrapped DNA, many tiny records, long contigs, protein FASTA,
  raw/gzip/BGZF transport, memory rows, and installed command-line comparators.

The benchmark binary generates deterministic synthetic FASTQ or FASTA in memory,
then measures the same parser and side-channel APIs used by downstream crates.
It reports best-of-N wall time to reduce noise from scheduler spikes.

For real input files, pass `--input PATH` or set `MICRORAPTOR_INPUT`. The file
path goes through `open_fastq_with_config`, so raw FASTQ, gzip FASTQ, and BGZF
FASTQ use the same auto-detection path as library callers.

For FASTA parser benchmarks, pass `--format fasta --mode parse`. FASTA rows use
`FastaReader` and the same raw/gzip/BGZF transport detection, but they are
parse-only: paired FASTQ validation and FASTQ pack rows do not apply.

Compression rows must be read literally. `flate2` and `libdeflate` are
third-party compression implementations; microraptor uses them as transport
backends. Benchmark rows named `libdeflate` measure microraptor parsing or
orchestration after selecting the upstream libdeflate engine through the
`libdeflater` Rust wrapper. They are not claims that microraptor implements a
new DEFLATE codec.

The FASTQ parser expects the common four-line record shape: name, sequence,
plus, quality. Multiline sequence or quality fields are not supported.

Paired input coverage is ordered-pair coverage: generated R1/R2 files with
matching order and interleaved files with adjacent mates. The stateful paired
reader handles different batch boundaries and validates normalized IDs, but it
does not synchronize reordered mates.

## Fast Commands

```bash
cargo test --all
cargo +nightly clippy --all-targets --all-features -- -D warnings
cargo +nightly clippy --all-targets --no-default-features -- -D warnings
cargo +nightly bench --all-features
scripts/bench.sh
scripts/prepare-real-benchmark-inputs.sh
scripts/benchmark-gauntlet.sh
scripts/render-benchmark-report.sh
scripts/check-benchmark-snapshots.sh
scripts/benchmark-common.sh
scripts/benchmark-rust-peers.sh
scripts/benchmark-fasta-peers.sh
scripts/benchmark-fasta-peer-size-sweep.sh
scripts/benchmark-fasta-gauntlet.sh
scripts/check-replication-host.sh --strict
scripts/check-pack-regression.sh
scripts/check-pack-instructions.sh
scripts/check-slab-autotune.sh
scripts/asm-pack.sh
scripts/release-gate.sh --allow-dirty
scripts/export-replication-kit.sh
```

For machine-readable output:

```bash
cargo run --release --bin microraptor-bench -- --records 500000 --iters 7 --json
cargo run --release --bin microraptor-bench -- --records 500000 --mode parse --json
cargo run --release --bin microraptor-bench -- --records 500000 --mode pack --json
cargo run --release --bin microraptor-bench -- --format fasta --mode parse --records 500000 --json
```

Synthetic `--mode pack` uses the trusted streaming pack path for `pack-seq-qual`.
It also reports `direct-pack-seq-qual` for the lower-memory single-pass scanner
and `reader-pack-seq-qual` for the safe parser-backed reference. The default
trusted path reuses the SIMD newline scanner, handles slab carry and CRLF
trimming, uses the fused base+quality pack kernel, and skips batch record
construction before packing. The fused kernel packs four classified bases per
compact LUT lookup and accumulates qualities during the same exact-length walk.

Use `scripts/check-pack-regression.sh` as a narrow guard for the default pack
path. It checks that `pack-seq-qual` and `reader-pack-seq-qual` checksums match
and fails when the trusted path is more than the configured tolerance slower
than the reader-backed reference. It also checks that the direct scanner emits
the same checksum. In CI (`CI=true`), timing failures are disabled by default so
the gate remains checksum/shape stable on noisy runners; set
`MICRORAPTOR_ENFORCE_TIMING=1` to make CI enforce the wall-clock threshold.

Use `scripts/check-bgzf-pack-regression.sh` as the matching guard for real BGZF
pack inputs. It creates both a small cyclic BGZF fixture and a large entropy
BGZF fixture that must cross the adaptive parallel threshold, requires the
serial/adaptive/libdeflate trusted-pack rows, checks their checksums against the
reader-backed reference, fails when the adaptive trusted BGZF path is more than
the configured tolerance slower than the reader-backed pack row, and checks that
the default BGZF pack row stays within tolerance of the explicit adaptive row.
The shell script only orchestrates fixtures; `microraptor-bench
--check-bgzf-pack-regression` owns row validation and tolerance checks. It uses
the same CI timing policy as the raw pack guard: timing checks are skipped under
`CI=true` unless `MICRORAPTOR_ENFORCE_TIMING=1`.

Use `scripts/asm-pack.sh` to emit optimized assembly for pack-path inspection.
Use `scripts/check-pack-instructions.sh` to compute pack-path
instructions/base from `perf stat`; set
`MICRORAPTOR_MAX_PACK_INSTRUCTIONS_PER_BASE` to turn it into a threshold gate.
Use `scripts/check-slab-autotune.sh` to compare the default pack path across
candidate slab sizes. Set `MICRORAPTOR_SLAB_INPUT=/path/to/file.fastq` to tune a
real workload instead of the synthetic fixture. Treat this as workload evidence,
not an automatic default change; slab winners are sensitive to input size,
compression, and cache state.

Set `MICRORAPTOR_GAUNTLET_CORPUS_INPUTS` to a space-separated list of real FASTQ,
gzip FASTQ, or BGZF FASTQ files to add optional corpus rows to the gauntlet
without checking datasets into the repository. These external corpus rows are
evidence capture, not a hard repository gate: successful rows are appended to the
JSONL output, and failed rows stay in the markdown report with their exit status
so malformed or unsupported local datasets do not hide the synthetic gauntlet
result. Single-file corpus inputs also get command-line comparator rows for
`seqkit`, `seqtk`, `samtools import -0`, and single-end `fastp` when those tools
are available on `PATH`.

Set `MICRORAPTOR_GAUNTLET_CORPUS_PAIRED_INPUTS` to a space-separated list of
`R1,R2,label` triples to add real paired-end rows. The label is optional but
recommended for stable report names. These paired corpus rows are also included
in external comparator commands when `seqkit`, `seqtk`, `samtools`, or `fastp`
are available on `PATH`.

Recommended corpus pass:

```bash
MICRORAPTOR_GAUNTLET_RESULT_DIR=target/bench-results-real \
MICRORAPTOR_GAUNTLET_CORPUS_INPUTS="/path/to/r1.fastq.gz /path/to/r2.fastq.gz /path/to/reads.fastq.bgz" \
MICRORAPTOR_GAUNTLET_CORPUS_PAIRED_INPUTS="/path/to/r1.fastq.gz,/path/to/r2.fastq.gz,real-r1-r2" \
scripts/benchmark-gauntlet.sh
```

Keep the corpus list outside git. For release evidence, use at least one real
short-read R1/R2 gzip pair and one BGZF file large enough to cross the adaptive
parallel threshold. If a local file fails parsing, keep the markdown report; the
failure is useful compatibility evidence but should not be counted as biological
throughput proof.

Prepare local Drosophila gzip/BGZF derivatives with:

```bash
PATH=~/miniconda3/envs/bench/bin:$PATH scripts/prepare-real-benchmark-inputs.sh
```

This writes `target/bench-real-inputs/drosophila_melanogaster/manifest.tsv`.
The combined BGZF derivative is intentionally built from the real R1 and R2
FASTQ files to exceed the adaptive BGZF threshold while keeping the source data
outside git.

The Rust parser-library peer harness reads one in-memory FASTQ byte buffer
through Microraptor, `seq_io`, `noodles-fastq`, and `bio`. It reports
`microraptor-stream` for the validated batch reader over a `Read` source and
`microraptor-slice-visitor` for the validated zero-copy visitor over a resident
byte slice. Its default consumer mode is `light`, which records record/base
accounting and a small shape checksum so parser framing is not hidden by hashing
every base. Set `MICRORAPTOR_RUST_PEER_CONSUMER=full` to hash every sequence
and quality byte. Set `MICRORAPTOR_RUST_PEER_DIAGNOSTICS=1` to include
additional visitor, trusted/no-validation, and raw record-ref rows.

Discover the broader local biological corpus with:

```bash
scripts/discover-local-benchmark-corpus.sh
```

The script reads `~/Projects/Benchmarks` by default and writes
`target/bench-corpus/local-corpus.tsv` plus
`target/bench-corpus/recommended-gauntlet.env`. The recommended env file
includes the bounded Drosophila ONT/PacBio single-end rows and the 1M paired
Illumina row when present. It also writes
`target/bench-corpus/independent-gauntlet.env` for bounded non-Drosophila E.
coli and yeast paired rows, and `target/bench-corpus/larger-gauntlet.env` for
explicit larger paired-end rows such as 5M. See `docs/benchmarks/CORPUS.md` for
the corpus ladder and stress run boundaries.

The checked local Drosophila snapshot was generated from the prepared benchmark
corpus with:

```bash
PATH=~/miniconda3/envs/bench/bin:$PATH \
MICRORAPTOR_GAUNTLET_RESULT_DIR=target/bench-results-drosophila-1m \
MICRORAPTOR_GAUNTLET_CORPUS_PAIRED_INPUTS=~/Projects/Benchmarks/datasets/drosophila_melanogaster/illumina_pe_r1.1m.fq,~/Projects/Benchmarks/datasets/drosophila_melanogaster/illumina_pe_r2.1m.fq,drosophila_illumina_1m \
scripts/benchmark-gauntlet.sh

scripts/render-benchmark-report.sh \
  target/bench-results-drosophila-1m/microraptor-gauntlet.jsonl \
  docs/benchmarks/drosophila-1m
```

The rendered artifact is
`docs/benchmarks/drosophila-1m/summary.md`. The same directory contains the
sanitized raw `microraptor-gauntlet.jsonl` used to render the tables. It
includes 2,000,000 FASTQ records, 72,000,000 sequenced bases, tool versions,
exact comparator commands, and SVG figures generated by repository scripts.

The checked compressed Drosophila snapshot was generated with:

```bash
PATH=~/miniconda3/envs/bench/bin:$PATH \
MICRORAPTOR_GAUNTLET_RESULT_DIR=target/bench-results-drosophila-compressed \
MICRORAPTOR_GAUNTLET_CORPUS_INPUTS=target/bench-real-inputs/drosophila_melanogaster/illumina_pe_r1_r2.1m-combined.fq.bgz \
MICRORAPTOR_GAUNTLET_CORPUS_PAIRED_INPUTS="target/bench-real-inputs/drosophila_melanogaster/illumina_pe_r1.1m.fq.gz,target/bench-real-inputs/drosophila_melanogaster/illumina_pe_r2.1m.fq.gz,drosophila_illumina_1m_gzip target/bench-real-inputs/drosophila_melanogaster/illumina_pe_r1.1m.fq.bgz,target/bench-real-inputs/drosophila_melanogaster/illumina_pe_r2.1m.fq.bgz,drosophila_illumina_1m_bgzf" \
scripts/benchmark-gauntlet.sh

cp target/bench-real-inputs/drosophila_melanogaster/manifest.tsv \
  docs/benchmarks/drosophila-compressed/input-manifest.tsv

scripts/render-benchmark-report.sh \
  target/bench-results-drosophila-compressed/microraptor-gauntlet.jsonl \
  docs/benchmarks/drosophila-compressed
```

The rendered artifact is
`docs/benchmarks/drosophila-compressed/summary.md`. The same directory contains
the sanitized raw `microraptor-gauntlet.jsonl` used to render the tables. It
includes a real paired gzip pass and a combined real BGZF input with
`54,594,251` compressed bytes, above the default `32 MiB` adaptive parallel
threshold.

The checked Drosophila read-type snapshot was generated with:

```bash
PATH=~/miniconda3/envs/bench/bin:$PATH \
scripts/discover-local-benchmark-corpus.sh

source target/bench-corpus/recommended-gauntlet.env

PATH=~/miniconda3/envs/bench/bin:$PATH \
MICRORAPTOR_GAUNTLET_RESULT_DIR=target/bench-results-drosophila-read-types \
scripts/benchmark-gauntlet.sh

cp target/bench-corpus/local-corpus.tsv \
  docs/benchmarks/drosophila-read-types/input-manifest.tsv

scripts/render-benchmark-report.sh \
  target/bench-results-drosophila-read-types/microraptor-gauntlet.jsonl \
  docs/benchmarks/drosophila-read-types
```

The rendered artifact is
`docs/benchmarks/drosophila-read-types/summary.md`. The same directory contains
the sanitized raw `microraptor-gauntlet.jsonl` used to render the tables. It
includes real Drosophila Illumina PE, PacBio CLR, and ONT FASTQ rows plus
installed `seqkit`, `seqtk`, `samtools`, and `fastp` comparator timings. The 5M
paired-end row is excluded from the recommended env because external workflow
comparators can dominate wall time; source
`target/bench-corpus/larger-gauntlet.env` explicitly when that cost is intended.

The checked independent-organism snapshot was generated with:

```bash
PATH=~/miniconda3/envs/bench/bin:$PATH \
scripts/discover-local-benchmark-corpus.sh

source target/bench-corpus/independent-gauntlet.env

PATH=~/miniconda3/envs/bench/bin:$PATH \
MICRORAPTOR_GAUNTLET_RESULT_DIR=target/bench-results-independent-organisms \
scripts/benchmark-gauntlet.sh

cp target/bench-corpus/local-corpus.tsv \
  docs/benchmarks/independent-organisms/input-manifest.tsv

scripts/render-benchmark-report.sh \
  target/bench-results-independent-organisms/microraptor-gauntlet.jsonl \
  docs/benchmarks/independent-organisms
```

The rendered artifact is
`docs/benchmarks/independent-organisms/summary.md`. The same directory contains
the sanitized raw `microraptor-gauntlet.jsonl` used to render the tables. It
includes real E. coli MG1655 and Saccharomyces cerevisiae BTT paired-end FASTQ
rows plus installed `seqkit`, `seqtk`, `samtools`, and `fastp` comparator
timings. This is independent organism coverage on the same workstation, not
independent-machine replication.

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
- `target/bench-results/microraptor-gauntlet-metadata.md`
- `target/bench-results/external-tools.tsv`

It also records whether optional external comparators such as `seqkit` or
`fastp` were installed and runnable. Metadata includes git commit and dirty
state, Rust and Cargo versions, release/all-feature build mode, CPU, logical CPU
count, RAM, kernel, filesystem class, available storage, comparator versions,
and gauntlet parameters.

Render the Microraptor auto-parse summary and SVG figure from the gauntlet JSONL:

```bash
scripts/render-benchmark-report.sh
```

By default this writes:

- `docs/benchmarks/latest/summary.md`
- `docs/benchmarks/latest/microraptor-gauntlet.jsonl`
- `docs/benchmarks/latest/metadata.md` when gauntlet metadata is available
- `docs/benchmarks/latest/external-tools.tsv` when external timings are available
- `docs/benchmarks/latest/figures/auto-bases-throughput.svg`

Use an explicit input/output pair for real-data snapshots:

```bash
scripts/render-benchmark-report.sh \
  target/bench-results-real/microraptor-gauntlet.jsonl \
  docs/benchmarks/real-snapshot
```

Verify that checked gauntlet and Rust peer summaries and SVG figures still
match their stored raw JSONL/TSV inputs with:

```bash
scripts/check-benchmark-snapshots.sh
```

Export a replication bundle for another machine with:

```bash
scripts/export-replication-kit.sh
```

See `docs/REPLICATION.md` for independent-machine evidence requirements.

Run Rust parser-library peer comparisons separately from the gauntlet:

```bash
MICRORAPTOR_RUST_PEER_ITERS=3 scripts/benchmark-rust-peers.sh

MICRORAPTOR_RUST_PEER_ITERS=3 \
MICRORAPTOR_RUST_PEER_OUT_DIR=docs/benchmarks/rust-peers-drosophila-r1 \
MICRORAPTOR_RUST_PEER_INPUT=~/Projects/Benchmarks/datasets/drosophila_melanogaster/illumina_pe_r1.1m.fq \
scripts/benchmark-rust-peers.sh
```

The script writes `rust-library-peers.tsv`, `summary.md`, `metadata.md`, and an
SVG throughput figure. It generates a temporary Cargo project under
`target/rust-peer-bench` so peer crates do not become normal crate dependencies.
The comparison is raw FASTQ parser-library evidence only; do not use it to make
claims about gzip, BGZF, command-line preprocessing, trimming, or filtering
behavior.

For a publication-style parser-framework scaling plot with time on the Y axis
and input size on the X axis:

```bash
MICRORAPTOR_SIZE_SWEEP_RECORDS="10000 50000 100000 500000 1000000 5000000" \
MICRORAPTOR_SIZE_SWEEP_COMPRESSIONS="raw gzip" \
MICRORAPTOR_SIZE_SWEEP_ITERS=5 \
scripts/benchmark-rust-peer-size-sweep.sh
```

This writes:

- `target/bench-results/rust-peer-size-sweep/rust-peer-size-sweep.tsv`
- `target/bench-results/rust-peer-size-sweep/summary.md`
- `target/bench-results/rust-peer-size-sweep/metadata.md`
- `target/bench-results/rust-peer-size-sweep/figures/rust-peer-size-sweep-time.svg`

The sweep reuses `scripts/benchmark-rust-peers.sh` for each input size and
compression mode, so every row still checks count/checksum parity across
`microraptor`, `seq_io`, `noodles-fastq`, and `bio`. Raw rows include the
resident `microraptor-slice-visitor`; gzip rows exclude it because that API is
intentionally raw resident-byte only. Treat the sweep as parser-framework
evidence; run the gauntlet separately for command-line workflow comparisons.

Run FASTA parser-library peer comparisons separately from the FASTQ gauntlet:

```bash
MICRORAPTOR_BENCH_THREADS=8 MICRORAPTOR_FASTA_PEER_ITERS=3 \
scripts/benchmark-fasta-peers.sh
```

The FASTA peer harness writes `fasta-library-peers.tsv`, `summary.md`,
`metadata.md`, and an SVG throughput figure. It compares robust multiline
resident parsing, strict two-line resident visitors, strict two-line streaming
visitors, strict resident/streaming light-accounting two-line counter paths,
`seq_io`, and `bio`. Rows name whether they use raw bytes, third-party flate2
gzip, or third-party libdeflate gzip via `libdeflater`. Treat strict `two-line`
rows as evidence for canonical `>header`/`sequence` FASTA, not arbitrary
multiline FASTA.

The FASTQ and FASTA peer harnesses share `scripts/benchmark-common.sh` for the
8-thread benchmark cap, sanitized checked artifacts, display paths, and generated
Cargo dependency specifications. Keep new benchmark families on that helper
instead of adding bespoke shell glue.

For a FASTA publication-style parser-framework scaling plot:

```bash
MICRORAPTOR_BENCH_THREADS=8 \
MICRORAPTOR_FASTA_SIZE_SWEEP_RECORDS="10000 50000 100000 500000 1000000" \
MICRORAPTOR_FASTA_SIZE_SWEEP_COMPRESSIONS="raw gzip" \
MICRORAPTOR_FASTA_SIZE_SWEEP_ITERS=3 \
scripts/benchmark-fasta-peer-size-sweep.sh
```

The checked local artifact lives at
[`docs/benchmarks/fasta-peer-size-sweep/summary.md`](docs/benchmarks/fasta-peer-size-sweep/summary.md).
It includes the raw TSV and
[`figures/fasta-peer-size-sweep-time.svg`](docs/benchmarks/fasta-peer-size-sweep/figures/fasta-peer-size-sweep-time.svg).
By default this script fails if the best row for any requested size/compression
does not start with `microraptor`; set
`MICRORAPTOR_FASTA_SIZE_SWEEP_REQUIRE_MICRORAPTOR_WINS=0` only for diagnostics.

For a broader FASTA robustness/performance gauntlet:

```bash
MICRORAPTOR_BENCH_THREADS=8 \
MICRORAPTOR_FASTA_GAUNTLET_RECORDS=10000 \
MICRORAPTOR_FASTA_GAUNTLET_READ_LEN=150 \
MICRORAPTOR_FASTA_GAUNTLET_ITERS=3 \
scripts/benchmark-fasta-gauntlet.sh
```

The checked local artifact lives at
[`docs/benchmarks/fasta-gauntlet/summary.md`](docs/benchmarks/fasta-gauntlet/summary.md).
It covers synthetic two-line DNA, wrapped DNA, many tiny records, long wrapped
contigs, protein FASTA, raw/gzip/BGZF transport, RSS smoke rows, optional local
corpus FASTA files, and installed command-line comparator timings.

`scripts/check-benchmark-snapshots.sh` verifies the checked FASTA size-sweep
summary/figure from the TSV and metadata, verifies the FASTA gauntlet summary
and required sidecar TSV files, and rejects legacy shell-wrapped external
commands across all checked `external-tools.tsv` artifacts.

Use the opt-in diagnostic mode to investigate microraptor internals without
changing the normal published peer table:

```bash
MICRORAPTOR_RUST_PEER_DIAGNOSTICS=1 scripts/benchmark-rust-peers.sh

MICRORAPTOR_RUST_PEER_CARGO='cargo' \
MICRORAPTOR_RUST_PEER_MICRORAPTOR_FEATURES=simd \
MICRORAPTOR_RUST_PEER_INPUT=~/Projects/Benchmarks/datasets/drosophila_melanogaster/illumina_pe_r1.1m.fq \
scripts/benchmark-rust-peers.sh
```

The diagnostic rows split validated microraptor parsing from no-validation and
direct `record_refs` iteration. On the local Drosophila R1 row, those variants
showed that validation and accessors were not the main loss; newline discovery
was. Stable default newline search now uses `memchr`; the `simd` feature uses
stable x86_64 `std::arch` acceleration when AVX2 is available and falls back to
scalar code elsewhere.

Build with `--features libdeflate` or `--all-features` to include explicit
`bgzf-libdeflate-*` rows for synthetic BGZF and `file-bgzf-libdeflate-*` rows
for real `.bgz` inputs. Ordinary gzip remains on the streaming flate2 path;
flate2 is configured for its pure-Rust backend by default.
`open_fastq_gzip_libdeflate` is an explicit buffered path for bounded gzip
inputs. Both flate2 and libdeflate are third-party compression backends. The
normal BGZF auto-open path uses `BgzfAutoReader` by default: small compressed
inputs stay serial and larger inputs switch to the bounded parallel reader only
past the built-in size threshold. When the `libdeflate` feature is enabled,
BGZF auto-open uses the upstream libdeflate inflate backend by default; use
`open_fastq_bgzf_flate2` or `open_fastq_bgzf_with_backend` when comparing or
forcing a backend. BGZF output can use libdeflate through `BgzfDeflateBackend`.
Use `open_fastq_bgzf_adaptive` when you need to pass an explicit BGZF worker or
backend configuration. `BgzfParallelConfig::with_parallel_min_compressed_bytes`
controls the adaptive serial/parallel threshold; `microraptor-bench` exposes the
same knob as `--bgzf-parallel-min-bytes`.

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
MICRORAPTOR_PROFILE_INPUT=reads.fastq.bgz MICRORAPTOR_PROFILE_BGZF_PARALLEL=1 scripts/profile-hotpath.sh
```

`profile-hotpath.sh` writes:

- `target/profiles/microraptor-hotpath.jsonl`
- `target/profiles/microraptor-parse.perf-stat.txt` when `perf` is permitted
- `target/profiles/microraptor-pack.perf-stat.txt` when `perf` is permitted

For BGZF inputs, `MICRORAPTOR_PROFILE_BGZF_PARALLEL=1` adds
`file-bgzf-direct-parallel-pack-seq-qual` to the pack benchmark so adaptive
BGZF overhead can be compared directly against a forced `BgzfParallelReader`.
BGZF adaptive and forced-parallel JSON rows include
`bgzf_job_queue_full` and `bgzf_result_queue_full` counters. These count bounded
queue backpressure events in the decompress/parse/pack path and are meant to
distinguish reader-starved, worker-starved, and consumer-starved profiles from
plain wall-clock noise.

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
In pack mode, `file-pack-seq-qual` uses the default adaptive BGZF auto-open path
for `.bgz` inputs. `bgzf-adaptive-pack-seq-qual` and
`file-bgzf-adaptive-pack-seq-qual` keep the adaptive comparison explicit. With
`libdeflate` enabled, the `*-libdeflate-serial-pack-seq-qual` and
`*-libdeflate-adaptive-pack-seq-qual` rows make the backend comparison explicit.

## CI Parity

The GitHub Actions workflow splits crate-readiness from the nightly performance
surface.

Stable crate surface:

```bash
cargo fmt --all -- --check
cargo clippy --lib -- -D warnings
cargo test --lib
cargo package
```

Nightly all-feature surface:

```bash
cargo +nightly clippy --all-targets --all-features -- -D warnings
cargo +nightly clippy --all-targets --no-default-features -- -D warnings
cargo +nightly test --all-features
cargo +nightly test --no-default-features
cargo +nightly fuzz build
```

`cargo fuzz build` only compiles the fuzz targets. It does not run long fuzzing
campaigns in CI.

## Publication Protocol

A public benchmark claim needs more than a local timing table. Before publishing
numbers, capture:

- Microraptor git commit and feature flags.
- Rust stable and nightly versions.
- CPU model, physical/logical core count, RAM, kernel, OS, and storage class.
- Input source, read length distribution, compression format, compressed bytes,
  records, and bases.
- Comparator versions and exact commands.
- Raw JSONL, rendered summary, and SVG figures generated by repository scripts.
- Checksum or count parity for every row used in a comparison.
- Separation between command-line workflow comparisons and Rust parser-library
  comparisons.

Minimum matrix for a paper-style result:

- Synthetic single-end raw/gzip/BGZF.
- Synthetic paired R1/R2 raw/gzip/BGZF.
- Synthetic interleaved raw/gzip/BGZF.
- At least one real short-read paired gzip dataset.
- At least one real BGZF dataset above the adaptive parallel threshold.
- Comparator rows for available command-line tools: `fastp`, `seqkit`, `seqtk`,
  `samtools import`, and `bgzip` where each tool meaningfully supports the task.
- Rust parser-library peer rows for `seq_io`, `noodles-fastq`, and `bio`.

Do not compare compressed `input_mib_s` against raw `input_mib_s` as a biological
throughput claim. Prefer `records_s` and `bases_s`, and state whether the row is
parse-only, pack-only, or parse-plus-pack.
