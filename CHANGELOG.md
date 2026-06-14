# Changelog

All notable changes to microraptor are documented here. The project is pre-1.0;
public APIs may still change between minor releases, but release notes should
call out any behavior, feature-flag, or benchmark-surface changes that affect
scientific users.

## Unreleased

### Added

- Stable default crate surface for raw, gzip, and BGZF FASTQ streaming.
- Ordered paired-end parsing for separate R1/R2 streams and adjacent
  interleaved records.
- Single-pass FASTQ visitor APIs: `FastqReader::visit_records` for streaming
  readers and `visit_fastq_bytes` for complete resident FASTQ byte buffers.
- Streaming multiline FASTA parsing via `FastaReader`, `FastaBatch`,
  `FastaRecord`, `FastaReader::visit_records`, `open_fasta`, and
  `open_fasta_with_config` over raw, gzip, and BGZF transport.
- Resident FASTA visitors via `visit_fasta_bytes`, plus strict two-line FASTA
  fast paths via `visit_two_line_fasta_bytes` and `visit_two_line_fasta_read`.
- FASTA shape detection, automatic resident fast-path dispatch, and strict
  two-line count/stat APIs via `detect_fasta_shape`, `visit_fasta_bytes_auto`,
  `count_two_line_fasta_bytes`, `count_two_line_fasta_read`, and `FastaStats`.
- Robust multiline FASTA stats APIs via `FastaReader::stats`,
  `count_fasta_read`, and `count_fasta_bytes`.
- `.fai`-style FASTA reference indexing via `FastaIndex`, `FastaIndexEntry`,
  `build_fasta_index`, and BGZF sequence-start virtual-offset annotation via
  `build_fasta_index_bgzf`.
- FASTA fixture generation for two-line/wrapped layouts and DNA/protein
  alphabets in `microraptor-fixture`.
- Configurable pair validation modes: full normalized ID validation, fast
  `/1`/`/2` validation, and trusted no-validation mode.
- Trusted four-line FASTQ pack paths for packed two-bit bases, ambiguity masks,
  Phred+33 summaries, and paired packed records.
- BGZF reader, writer, adaptive serial/parallel reader, parallel whole-buffer
  helpers, virtual-offset index construction, and seek reader.
- Optional `libdeflate` BGZF inflate/deflate backends and explicit buffered
  FASTQ/FASTA gzip openers.
- FASTA fuzz target covering streaming reader, resident visitors, robust stats,
  shape detection, and index construction.
- Benchmark gauntlet with synthetic raw/gzip/BGZF fixtures, optional local
  corpus rows, external command-line comparator rows, and rendered Markdown/SVG
  reports.
- Checked benchmark snapshot verifier that re-renders stored gauntlet JSONL and
  Rust/FASTA peer TSV snapshots, then diffs generated summaries/figures against
  checked artifacts; FASTA gauntlet and FASTA size-sweep artifacts are now part
  of the verifier.
- Local release-gate script for clean-tree package/readiness checks, with
  optional nightly and benchmark regeneration surfaces.
- Replication protocol and replication-kit exporter for independent-machine
  benchmark review.
- Rust parser-library peer benchmark script for `seq_io`, `noodles-fastq`, and
  `bio`.
- FASTA parse-only benchmark mode in `microraptor-bench` via
  `--format fasta --mode parse`.
- FASTA Rust peer benchmark scripts and a checked raw/gzip size-sweep artifact
  under `docs/benchmarks/fasta-peer-size-sweep/`.
- FASTA shape/transport gauntlet under `scripts/benchmark-fasta-gauntlet.sh`
  with checked artifacts under `docs/benchmarks/fasta-gauntlet/`.
- Shared benchmark harness helpers in `scripts/benchmark-common.sh`, shared
  fixture generation in `benchutil`, and a split `microraptor-bench` report
  module.
- Local biological corpus discovery for `~/Projects/Benchmarks`, including
  Drosophila Illumina PE, PacBio CLR, and ONT read-type rows.
- GitHub issue and pull request templates for parser bugs, benchmark claims,
  feature requests, and review evidence.
- Security policy for parser, pack, and BGZF vulnerability reports.
- Crate-level rustdoc examples and warning-denied missing-doc coverage for the
  public reader, opener, pairing, error, pack, and BGZF surfaces.

### Fixed

- FASTQ slab framing now preserves a trailing partial line after complete
  records at slab boundaries.
- Fast `/1`/`/2` pair validation rejects wrong mate suffix combinations.
- Interleaved paired reads now honor the configured pair-validation mode.
- Stable default FASTQ newline discovery now uses `memchr`, closing most of the
  parser-only gap to Rust peer libraries while keeping nightly SIMD opt-in.
- FASTQ streaming batch framing now pre-reserves record side tables and avoids
  per-record fallible `u32` range conversions after a slab-level bounds check.
- Strict two-line FASTA stream visitor and stream counter now share one scanner
  implementation.

### Benchmark Evidence

- Synthetic fixture snapshot:
  `docs/benchmarks/latest/summary.md`.
- Real Drosophila 1M paired raw snapshot:
  `docs/benchmarks/drosophila-1m/summary.md`.
- Real Drosophila-derived gzip/BGZF snapshot:
  `docs/benchmarks/drosophila-compressed/summary.md`.
- Real Drosophila read-type snapshot covering Illumina PE, PacBio CLR, and ONT:
  `docs/benchmarks/drosophila-read-types/summary.md`.
- Independent-organism paired FASTQ snapshot covering E. coli MG1655 and yeast
  BTT:
  `docs/benchmarks/independent-organisms/summary.md`.
- Rust parser-library peer snapshots:
  `docs/benchmarks/rust-peers/summary.md` and
  `docs/benchmarks/rust-peers-drosophila-r1/summary.md`. These snapshots now
  separate the validated streaming batch row from the validated resident-slice
  visitor row and use light parser accounting by default; set
  `MICRORAPTOR_RUST_PEER_CONSUMER=full` to hash every sequence and quality byte.
- FASTA parser-framework size sweep:
  `docs/benchmarks/fasta-peer-size-sweep/summary.md`, generated with
  `MICRORAPTOR_BENCH_THREADS=8` and a required microraptor-winner guard for
  each requested raw/gzip size.
- FASTA shape/transport gauntlet:
  `docs/benchmarks/fasta-gauntlet/summary.md`, covering two-line DNA, wrapped
  DNA, many tiny records, long wrapped contigs, protein FASTA, raw/gzip/BGZF
  transport, RSS smoke rows, and installed command-line comparator timings.

### Claim Boundary

- Current benchmark artifacts are local machine evidence. They support the
  existence and reproducibility of the benchmark machinery, but they should not
  be presented as universal speed claims without regenerating from the release
  commit and recording hardware, toolchain, comparator versions, and exact
  commands.
