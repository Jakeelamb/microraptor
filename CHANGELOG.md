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
- Configurable pair validation modes: full normalized ID validation, fast
  `/1`/`/2` validation, and trusted no-validation mode.
- Trusted four-line FASTQ pack paths for packed two-bit bases, ambiguity masks,
  Phred+33 summaries, and paired packed records.
- BGZF reader, writer, adaptive serial/parallel reader, parallel whole-buffer
  helpers, virtual-offset index construction, and seek reader.
- Optional `libdeflate` BGZF inflate/deflate backends and explicit buffered
  gzip opener.
- Benchmark gauntlet with synthetic raw/gzip/BGZF fixtures, optional local
  corpus rows, external command-line comparator rows, and rendered Markdown/SVG
  reports.
- Checked benchmark snapshot verifier that re-renders stored gauntlet JSONL and
  Rust peer TSV snapshots, then diffs generated summaries/figures against
  checked artifacts.
- Local release-gate script for clean-tree package/readiness checks, with
  optional nightly and benchmark regeneration surfaces.
- Replication protocol and replication-kit exporter for independent-machine
  benchmark review.
- Rust parser-library peer benchmark script for `seq_io`, `noodles-fastq`, and
  `bio`.
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
  `docs/benchmarks/rust-peers-drosophila-r1/summary.md`.

### Claim Boundary

- Current benchmark artifacts are local machine evidence. They support the
  existence and reproducibility of the benchmark machinery, but they should not
  be presented as universal speed claims without regenerating from the release
  commit and recording hardware, toolchain, comparator versions, and exact
  commands.
