# Replication Protocol

Microraptor benchmark artifacts are useful only when their scope is explicit.
The checked snapshots in this repository are local-machine evidence. They
demonstrate that the benchmark machinery works, that competitor rows can be
captured, and that raw artifacts reproduce checked summaries and figures. They
do not establish universal performance.

## Release-Commit Regeneration

From a clean release commit, run:

```bash
scripts/check-replication-host.sh --strict
scripts/release-gate.sh --nightly --bench
```

The preflight records host metadata and verifies the tools needed for release
and benchmark regeneration. The release gate verifies formatting,
warning-denied library and rustdoc surfaces, tests, release documentation,
script syntax, checked benchmark snapshots, nightly feature modes, fuzz target
compilation, synthetic benchmark regeneration, Rust peer benchmark regeneration,
and `cargo package`.

For real-data snapshots, use the local corpus scripts:

```bash
PATH=~/miniconda3/envs/bench/bin:$PATH \
scripts/discover-local-benchmark-corpus.sh

source target/bench-corpus/recommended-gauntlet.env

PATH=~/miniconda3/envs/bench/bin:$PATH \
MICRORAPTOR_GAUNTLET_RESULT_DIR=target/bench-results-real \
scripts/benchmark-gauntlet.sh

scripts/render-benchmark-report.sh \
  target/bench-results-real/microraptor-gauntlet.jsonl \
  docs/benchmarks/real-snapshot
```

If the host has the Trex-local E. coli and yeast corpus rows, also run:

```bash
source target/bench-corpus/independent-gauntlet.env

PATH=~/miniconda3/envs/bench/bin:$PATH \
MICRORAPTOR_GAUNTLET_RESULT_DIR=target/bench-results-independent-organisms \
scripts/benchmark-gauntlet.sh

scripts/render-benchmark-report.sh \
  target/bench-results-independent-organisms/microraptor-gauntlet.jsonl \
  docs/benchmarks/independent-organisms
```

Finish with:

```bash
scripts/check-benchmark-snapshots.sh
```

## Replication Kit

Create a portable evidence bundle with:

```bash
scripts/export-replication-kit.sh
```

The kit is written to `target/replication-kit` by default. It includes checked
benchmark summaries, raw JSONL/TSV artifacts, generated figures, release docs,
benchmark scripts, a replication README, and SHA-256 hashes for checked
artifacts. It does not include FASTQ datasets.

Use the kit to send a fixed evidence surface to another machine, then compare
newly generated raw JSONL/TSV and summaries against the checked artifacts and
claim boundaries.

Before running benchmarks on the local benchmark host, run:

```bash
scripts/check-replication-host.sh --strict
```

The script intentionally does not install dependencies. A missing core tool
such as `cargo` or `rustc` means that host is not yet usable for release
replication. Missing comparator tools mean the host can still exercise the Rust
crate, but cannot produce the full command-line competitor matrix.

## Independent-Machine Evidence

An independent-machine replication result should record:

- `scripts/check-replication-host.sh --strict` output;
- microraptor git commit and dirty/clean state;
- CPU model, logical CPU count, RAM, kernel, OS, filesystem, and storage class;
- Rust stable and nightly toolchains;
- feature flags and build profile;
- comparator versions for `seqkit`, `seqtk`, `samtools`, `bgzip`, and `fastp`;
- exact commands and environment variables;
- raw JSONL/TSV benchmark outputs;
- rendered Markdown summaries and SVG figures;
- checksum, record-count, and base-count parity for compared rows.

Acceptable public wording after one independent-machine run is still scoped:
"replicated on two machines for these datasets and commands." It is not a
claim of universal superiority across all FASTQ workloads or workflow tools.

## Claim Boundaries

Do not use replication artifacts to claim:

- support for multiline FASTQ;
- synchronization of reordered paired-end files;
- workflow equivalence between parser-only rows and tools that trim, filter,
  align, import, or emit QC reports;
- representative performance for every sequencing platform or machine.
