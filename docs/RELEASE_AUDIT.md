# Release Audit

This audit tracks whether microraptor is ready to present to the scientific
community. It is deliberately stricter than "tests pass": each release-facing
claim needs current evidence.

## Status

Current state: not final-release complete.

Microraptor has a publishable crate skeleton, a documented framework, local
benchmark snapshots, competitor comparisons, and generated figures. Remaining
work before a public performance announcement is mostly replication and final
release hygiene: regenerate all artifacts from the release commit, replicate the
real-data matrix on another machine, and rerun the final clean package gate from
the release commit.

## Requirement Matrix

| Requirement | Evidence | Status |
| --- | --- | --- |
| Crate metadata is suitable for crates.io | `Cargo.toml` has license, readme, repository, homepage, documentation, keywords, categories, and `rust-version`; `LICENSE-MIT` and `LICENSE-APACHE` are present | Satisfied |
| Default build is stable Rust | `rust-toolchain.toml` uses stable; default features are `bgzf` and `gzip`; nightly SIMD is behind `simd` | Satisfied |
| Public API is documented | `#![warn(missing_docs)]`; `RUSTFLAGS="-D warnings" cargo check --lib`; `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps` | Satisfied |
| Public API release surface is classified | `docs/API_SURFACE.md` tiers primary FASTQ readers/openers, advanced pack side channels, BGZF transport/indexing, and hidden bench helpers | Satisfied |
| README explains capabilities and limitations | `README.md` states scope, features, benchmark snapshots, limitations, examples, and claim boundary | Satisfied |
| Framework is ruthlessly analyzed | `docs/FRAMEWORK.md` defines the model, competitor categories, local evidence, and explicit "do claim" / "do not claim" boundaries | Satisfied |
| Independent replication protocol exists | `docs/REPLICATION.md`, `scripts/check-replication-host.sh`, and `scripts/export-replication-kit.sh` define host preflight, release-commit regeneration, replication-kit contents, independent-machine evidence fields, and claim boundaries | Satisfied |
| Benchmark protocol is reproducible | `BENCHMARKING.md`, `scripts/benchmark-gauntlet.sh`, `scripts/benchmark-common.sh`, `scripts/render-benchmark-report.sh`, `scripts/check-benchmark-snapshots.sh`, `scripts/benchmark-rust-peers.sh`, `scripts/benchmark-fasta-peers.sh`, `scripts/benchmark-fasta-peer-size-sweep.sh`, `scripts/benchmark-fasta-gauntlet.sh`, `scripts/prepare-real-benchmark-inputs.sh`, and `scripts/discover-local-benchmark-corpus.sh`; gauntlet metadata records commit state, toolchains, all-feature release mode, machine class, filesystem, comparator versions, and parameters | Satisfied |
| Benchmark figures and raw measurement artifacts are generated from scripts | `docs/benchmarks/*/figures/*.svg`, `docs/benchmarks/*/microraptor-gauntlet.jsonl`, `docs/benchmarks/*/microraptor-fasta-gauntlet.jsonl`, `docs/benchmarks/*/external-tools.tsv`, `docs/benchmarks/*/rust-library-peers.tsv`, and `docs/benchmarks/*/fasta-peer-size-sweep.tsv` rendered or copied by repository scripts with local home paths sanitized; `scripts/check-benchmark-snapshots.sh` re-renders checked gauntlet, Rust peer, FASTA peer, and FASTA size-sweep snapshots, checks FASTA gauntlet sidecars, diffs summaries/figures, and rejects legacy shell-wrapped external command records | Satisfied |
| Command-line competitor comparisons exist | `docs/benchmarks/drosophila-1m`, `docs/benchmarks/drosophila-compressed`, `docs/benchmarks/drosophila-read-types`, and `docs/benchmarks/independent-organisms` include `seqkit`, `seqtk`, `samtools`, and `fastp` rows where available | Satisfied |
| Rust parser-library peer comparisons exist | `docs/benchmarks/rust-peers` and `docs/benchmarks/rust-peers-drosophila-r1` compare `seq_io`, `noodles-fastq`, `bio`, and microraptor | Satisfied |
| Real biological datasets are covered | Local Drosophila Illumina PE, PacBio CLR, ONT, E. coli MG1655 paired-end, and yeast BTT paired-end rows are discovered and summarized; gzip/BGZF derivatives are benchmarked | Satisfied |
| Broad public performance claims are justified | Needs regeneration from release commit and independent-machine replication | Not satisfied |
| Citation metadata exists | `CITATION.cff` | Satisfied |
| Release notes exist | `CHANGELOG.md` | Satisfied |
| Contribution expectations are explicit | `CONTRIBUTING.md` | Satisfied |
| GitHub contribution workflow is structured | `.github/ISSUE_TEMPLATE/*`, `.github/pull_request_template.md`, and `SECURITY.md` request reproducible parser, benchmark, API, and vulnerability evidence | Satisfied |
| CI covers release-facing documentation gates | `.github/workflows/ci.yml` checks warning-denied library/rustdoc surfaces, release docs, benchmark script syntax, package dry-run, nightly feature modes, and fuzz target compilation | Satisfied |
| Final crates.io package verifies | `cargo package --allow-dirty` currently verifies; rerun on clean release tree before publishing | Partially satisfied |

## Current Local Verification

The current dirty development tree has passed:

```bash
scripts/release-gate.sh --allow-dirty --nightly
```

This covered formatting, warning-denied stable library check, warning-denied
rustdoc, `cargo test --all`, release artifact existence checks, benchmark script
syntax checks, checked benchmark snapshot verification, nightly all-feature
clippy, nightly no-default-feature clippy, nightly all-feature tests, nightly
no-default-feature tests, fuzz target compilation, and `cargo package
--allow-dirty`.

Observed test coverage from that run:

- Stable/default: 112 library tests, 5 benchmark-binary tests, 2 doctests.
- Nightly/all-features: 119 library tests, 5 benchmark-binary tests, 2 doctests.
- Nightly/no-default-features: 89 library tests, 5 benchmark-binary tests, 2
  doctests.
- Snapshot verifier: 5 gauntlet snapshots, 2 Rust peer snapshots, 0 checked
  standalone FASTA peer snapshots, 1 FASTA size-sweep snapshot, and 1 FASTA
  gauntlet snapshot.
- Package dry run: 160 files, 998.1 KiB unpacked, approximately 184.3 KiB
  compressed.

This is strong development evidence, but it is not a clean release-commit gate
because the current tree is intentionally dirty while release work is being
assembled.

## Required Release Gate

Run this from the release commit:

```bash
scripts/check-replication-host.sh --strict
scripts/release-gate.sh --nightly --bench
```

Expanded gate:

```bash
cargo fmt --all -- --check
cargo +stable test --lib
cargo +nightly test --all-features
cargo +nightly test --no-default-features
cargo +nightly clippy --all-targets --all-features -- -D warnings
cargo +nightly clippy --all-targets --no-default-features -- -D warnings
RUSTFLAGS="-D warnings" cargo check --lib
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps
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

For publishable benchmark figures, also regenerate the Drosophila real-data
snapshots from the release commit and verify that the checked raw JSONL/TSV
outputs match the rendered summaries.

## Current Claim Boundary

Safe claims:

- Microraptor is a Rust FASTQ streaming and packing core for raw, gzip, and BGZF
  inputs.
- It exposes borrowed FASTQ batches, ordered paired-end validation, BGZF
  transport helpers, and optional packed base/quality side channels.
- The repository contains reproducible benchmark scripts, local biological
  snapshots across Drosophila, E. coli, and yeast, command-line comparator rows,
  Rust parser-library peer rows, and generated figures.

Unsafe claims without more evidence:

- Faster than `fastp`, `seqkit`, `seqtk`, `samtools`, `seq_io`, noodles, or
  `bio` across workloads.
- Suitable for multiline FASTQ.
- Suitable for synchronizing reordered paired-end files.
- Representative of every sequencing platform or every machine.
