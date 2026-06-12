# Publishing Checklist

This repository should not be presented to the scientific community until the
items below are true and backed by current artifacts.

## Crate Readiness

- `Cargo.toml` has license, repository, homepage, documentation, README,
  keywords, categories, and `rust-version`.
- `LICENSE-MIT` and `LICENSE-APACHE` are present for the dual-license manifest.
- `CITATION.cff`, `CHANGELOG.md`, and `CONTRIBUTING.md` are present and current.
- Default features build on stable Rust.
- Nightly-only acceleration is behind explicit feature flags.
- Public examples compile in doc tests or unit tests.
- Public reader, opener, pairing, error, pack, and BGZF surfaces have rustdoc
  and pass the crate `missing_docs` lint with warnings denied.
- `docs/API_SURFACE.md` classifies the intended `0.1.x` public surface and
  identifies advanced APIs that should not be marketed as beginner entry points.
- `cargo package` succeeds on a clean tree before release.
- The README states limitations as clearly as capabilities.

## GitHub Readiness

- CI checks stable default library behavior.
- CI checks warning-denied library and rustdoc surfaces.
- CI checks nightly all-feature behavior, including no-default-features.
- Fuzz targets compile in CI.
- Issue templates distinguish parser/BGZF bugs, benchmark claims, and feature
  requests.
- The pull request template asks for checks, benchmark evidence, and public API
  documentation updates.
- `SECURITY.md` defines vulnerability-reporting scope for parser, pack, and
  BGZF behavior.
- Release notes list performance-affecting feature flags.
- Generated benchmark artifacts are reproducible from scripts.
- Release notes list current benchmark artifacts and claim boundaries.

## Scientific Readiness

- `docs/FRAMEWORK.md` defines the framework and claim boundaries.
- `docs/API_SURFACE.md` classifies the public API tiers and release-surface
  decision.
- `docs/REPLICATION.md` defines release-commit regeneration and
  independent-machine evidence requirements.
- `scripts/check-replication-host.sh --strict` passes on every host used for
  release or benchmark replication.
- `docs/RELEASE_AUDIT.md` maps release requirements to current evidence.
- `BENCHMARKING.md` gives an exact benchmark protocol and interpretation rules.
- Benchmark output includes raw JSONL and a rendered summary/figure.
- Comparator matrix includes command-line tools and Rust library peers.
- At least one real dataset pass is captured outside git and summarized with
  enough metadata to rerun.

## Benchmark Evidence Gate

Fast local release gate:

```bash
scripts/release-gate.sh
```

Use `scripts/release-gate.sh --allow-dirty` only for development verification
before the release commit exists. Use `--nightly` to include the full nightly
feature/fuzz surface and `--bench` to regenerate synthetic benchmark snapshots.

Minimum release evidence, expanded:

```bash
cargo +stable test --lib
cargo +nightly test --all-features
cargo +nightly test --no-default-features
cargo +nightly clippy --all-targets --all-features -- -D warnings
cargo +nightly clippy --all-targets --no-default-features -- -D warnings
RUSTFLAGS="-D warnings" cargo check --lib
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps
test -s CITATION.cff
test -s CHANGELOG.md
test -s CONTRIBUTING.md
test -s SECURITY.md
test -s .github/pull_request_template.md
test -s .github/ISSUE_TEMPLATE/parser_bug.yml
test -s .github/ISSUE_TEMPLATE/benchmark_claim.yml
test -s .github/ISSUE_TEMPLATE/feature_request.yml
test -s docs/API_SURFACE.md
test -s docs/REPLICATION.md
test -s docs/RELEASE_AUDIT.md
scripts/benchmark-gauntlet.sh
scripts/render-benchmark-report.sh
scripts/check-benchmark-snapshots.sh
scripts/benchmark-rust-peers.sh
scripts/check-replication-host.sh --strict
PATH=~/miniconda3/envs/bench/bin:$PATH scripts/prepare-real-benchmark-inputs.sh
PATH=~/miniconda3/envs/bench/bin:$PATH scripts/discover-local-benchmark-corpus.sh
scripts/export-replication-kit.sh
cargo package
```

Recommended real-data evidence:

```bash
source target/bench-corpus/recommended-gauntlet.env
PATH=~/miniconda3/envs/bench/bin:$PATH \
MICRORAPTOR_GAUNTLET_RESULT_DIR=target/bench-results-real \
scripts/benchmark-gauntlet.sh
scripts/render-benchmark-report.sh target/bench-results-real/microraptor-gauntlet.jsonl docs/benchmarks/real-snapshot
```

For publishable figures, record:

- CPU model and core count.
- RAM.
- Storage device or filesystem class.
- Kernel and OS.
- Rust toolchains.
- Microraptor git commit.
- Comparator versions.
- Exact commands and environment variables.

## Current Gaps

- The current Drosophila 1M local snapshot includes installed `seqkit`, `seqtk`,
  `samtools`, and `fastp` comparator rows, but it is still one machine-local
  dataset pass.
- The compressed Drosophila snapshot covers a real paired gzip pass and a real
  BGZF-derived input above the adaptive threshold, but it is still derived from
  the same Drosophila source dataset.
- The read-type Drosophila snapshot covers Illumina PE, PacBio CLR, and ONT
  FASTQ parser/comparator rows, but those long-read rows are parser
  compatibility evidence rather than assembly-quality evidence.
- The independent-organism snapshot covers E. coli MG1655 and yeast BTT paired
  FASTQ rows from the local Trex-governed corpus, but it is still same-machine
  evidence.
- External comparator rows are wall-clock command timings; public tables should
  still explain task differences between parser-only work and workflow tools.
- Existing target benchmark artifacts are useful for development; publishable
  figures should be regenerated from the release commit and paired with the
  checked rendered summary, raw JSONL, hardware metadata, and comparator
  versions.
- The real-data matrix now includes independent organisms, but it should still
  be replicated on at least one independent machine before broad public
  performance claims.
- Rust parser-library peer evidence now exists for synthetic raw FASTQ and
  Drosophila R1 raw FASTQ, but it should be expanded to compressed inputs only
  if the peer libraries are wrapped in equivalent decompression paths.
- Pack and BGZF public APIs now pass the crate `missing_docs` lint and are
  classified in `docs/API_SURFACE.md`; before `1.0`, revisit whether raw range
  types and whole-buffer BGZF helpers should remain as-is.
