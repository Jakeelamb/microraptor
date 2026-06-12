# Benchmark Corpus

Microraptor benchmark claims should separate repository fixtures from local
biological evidence. The repository fixtures are deterministic and cheap enough
for frequent regression checks. The local biological corpus is used for release
evidence and figure generation, but the FASTQ files themselves stay outside git.

## Local Workspace

On this workstation, `~/Projects/Benchmarks` currently provides the main real
FASTQ corpus:

| Dataset | Source | Size | Files | Role |
| --- | --- | --- | --- | --- |
| Drosophila Illumina PE 1M | ENA `DRR001444` | 1M read pairs / 72 Mbp | `illumina_pe_r1.1m.fq`, `illumina_pe_r2.1m.fq` | routine raw/gzip/BGZF release evidence |
| Drosophila Illumina PE 5M | ENA `DRR001444` | 5M read pairs / count verified by file presence, base count not in prepared manifest | `illumina_pe_r1.5m.fq`, `illumina_pe_r2.5m.fq` | explicit larger raw paired-end pass |
| Drosophila Illumina PE 25M | ENA `DRR001444` | 25M read pairs / 1.8 Gbp | `illumina_pe_r1.25m.fq`, `illumina_pe_r2.25m.fq` | deliberate stress runs only |
| Drosophila Illumina PE 50M | ENA `DRR001444` | 50M read pairs / 10 Gbp | `illumina_pe_r1.50m.fq`, `illumina_pe_r2.50m.fq` | deliberate stress runs only |
| Drosophila PacBio CLR 50k | SRA `SRR1204085` | 50k reads / 88.1 Mbp | `pacbio_clr_subreads.50k.fq` | long-read parser compatibility evidence |
| Drosophila ONT 50k | SRA `SRR18021217` | 50k reads / 173.2 Mbp | `ont.50k.fq` | long-read parser compatibility evidence |
| E. coli MG1655 Illumina PE 100k | ENA `SRR001666` plus RefSeq MG1655 reference | 100k read pairs | Trex-local `r1.100000.fq`, `r2.100000.fq` | non-Drosophila bacterial paired-end evidence |
| Yeast BTT Illumina PE 10k | ENA `ERR1308583` plus S288C reference | 10k read pairs | Trex-local `r1.10000.fq`, `r2.10000.fq` | non-Drosophila diploid eukaryotic paired-end evidence |

Counts and checksums for the Drosophila prepared files are in
`~/Projects/Benchmarks/manifests/drosophila_prepared.tsv`. Trex-local bacterial
and yeast rows are described in
`~/Projects/Benchmarks/manifests/local_datasets.tsv`.

Run discovery from the microraptor repository:

```bash
PATH=~/miniconda3/envs/bench/bin:$PATH scripts/check-replication-host.sh --strict
scripts/discover-local-benchmark-corpus.sh
```

This writes:

- `target/bench-corpus/local-corpus.tsv`: discovered FASTQ rows, roles, and
  sanitized paths.
- `target/bench-corpus/recommended-gauntlet.env`: shell variables for the
  recommended routine real-data gauntlet.
- `target/bench-corpus/independent-gauntlet.env`: shell variables for bounded
  E. coli and yeast paired-end rows used as independent-organism evidence.
- `target/bench-corpus/larger-gauntlet.env`: optional shell variables for larger
  paired-end rows that should be run deliberately.

## Recommended Release Pass

Use the benchmark conda environment for `seqkit`, `seqtk`, `samtools`, `bgzip`,
and `fastp`:

```bash
PATH=~/miniconda3/envs/bench/bin:$PATH \
scripts/discover-local-benchmark-corpus.sh

source target/bench-corpus/recommended-gauntlet.env

PATH=~/miniconda3/envs/bench/bin:$PATH \
MICRORAPTOR_GAUNTLET_RESULT_DIR=target/bench-results-drosophila-read-types \
scripts/benchmark-gauntlet.sh

scripts/render-benchmark-report.sh \
  target/bench-results-drosophila-read-types/microraptor-gauntlet.jsonl \
  docs/benchmarks/drosophila-read-types
```

This pass should include:

- Synthetic raw/gzip/BGZF fixtures for regression continuity.
- Real Drosophila Illumina paired-end raw FASTQ at 1M read pairs when present.
- Real Drosophila PacBio CLR and ONT single-end FASTQ rows.
- External command-line rows for `seqkit`, `seqtk`, `samtools`, and `fastp`
  where those tools are on `PATH`.

## Non-Drosophila Organism Pass

Use the independent env to add non-Drosophila local evidence from the Trex
governed corpus:

```bash
PATH=~/miniconda3/envs/bench/bin:$PATH \
scripts/discover-local-benchmark-corpus.sh

source target/bench-corpus/independent-gauntlet.env

PATH=~/miniconda3/envs/bench/bin:$PATH \
MICRORAPTOR_GAUNTLET_RESULT_DIR=target/bench-results-independent-organisms \
scripts/benchmark-gauntlet.sh

scripts/render-benchmark-report.sh \
  target/bench-results-independent-organisms/microraptor-gauntlet.jsonl \
  docs/benchmarks/independent-organisms
```

This pass should include E. coli MG1655 and Saccharomyces cerevisiae BTT paired
FASTQ rows when those Trex-local files are present. Treat it as cross-organism
evidence on this workstation, not as a substitute for a separate replication
run.

## Larger And Stress Passes

The 5M Drosophila paired-end row is useful larger-input evidence, but external
workflow comparators can dominate wall time. Source the larger env explicitly
after the recommended env when that cost is intended:

```bash
source target/bench-corpus/recommended-gauntlet.env
source target/bench-corpus/larger-gauntlet.env
```

The 25M and 50M Drosophila paired-end files are useful stress inputs, but they
are not the default release pass. They are large enough to dominate local disk
cache state and wall time, so run them explicitly and record that they were
stress runs.

## Claim Boundary

PacBio CLR and ONT rows are parser compatibility and throughput evidence only.
They do not establish assembly quality, read correction quality, or suitability
for every long-read workflow. Public claims still need exact command lines,
metadata, comparator versions, and generated figures from repository scripts.
