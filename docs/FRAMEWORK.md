# Microraptor Framework

Microraptor is a FASTQ streaming core, not a general sequencing workflow tool.
Its job is to turn raw, gzip, and BGZF FASTQ byte streams into validated,
borrowed record batches and optional packed base/quality side channels with the
least possible allocation and the clearest possible failure modes.

## Core Model

The framework is built around one invariant: decompression is a transport layer.
The parser consumes decompressed byte slabs regardless of whether the original
input was raw FASTQ, ordinary gzip, or BGZF. Public openers may choose different
transport backends, but downstream parsing and packing see the same byte model.

The current record model is intentionally narrow:

- FASTQ records are four physical lines: name, sequence, plus, quality.
- Records borrow from reusable slabs through `RecordRef` ranges.
- Batch lifetime is tied to the reader call, so the reader can reuse buffers.
- Quality and base packing can run as a side channel without forcing owned
  per-record allocation.
- Paired-end support assumes ordered mates, either in separate R1/R2 files or
  adjacent interleaved records.

That shape is the main advantage and the main constraint. Microraptor should be
excellent when a downstream tool needs a fast, low-allocation stream of ordinary
short-read FASTQ records. It should not pretend to replace tools that perform
adapter trimming, read filtering, quality-control reporting, alignment, or
format conversion.

## Execution Boundaries

### Stable Crate Surface

The default feature set is intended to build on stable Rust:

- `gzip`: gzip magic detection and streaming `flate2` input.
- `bgzf`: BGZF reader/writer, detection, virtual-offset indexing, and adaptive
  serial/parallel BGZF reading.

Nightly SIMD remains available behind the `simd` feature. This keeps crates.io
and downstream library adoption practical without removing the higher-throughput
research path.

Stable default newline discovery uses `memchr`; the older scalar byte loop was
the main reason microraptor trailed Rust parser peers on raw in-memory FASTQ
rows. Nightly SIMD is still useful: in the local Drosophila R1 peer harness, the
`simd` feature beat the checked `seq_io` row, while stable default narrowed the
gap to parser-library peers without requiring nightly.

### Compression Backends

Ordinary gzip auto-open uses a streaming `flate2::read::MultiGzDecoder`.
`open_fastq_gzip_libdeflate` is explicit because it buffers the decompressed
input and is only appropriate for bounded inputs.

BGZF is block-aware. The default BGZF auto reader uses serial reading for small
compressed inputs and bounded parallel reading above the adaptive threshold. With
`libdeflate`, BGZF inflate/deflate can use libdeflate-backed paths.

### Pairing Model

Microraptor validates ordered mates; it does not synchronize reordered files.
`PairValidation::Full` normalizes `/1` and `/2` suffixes before comparison.
`PairValidation::FastSlash` is a fast path for ordered `/1` and `/2` names with
fallback to full validation when the slash shape is not present.
`PairValidation::None` is for trusted internal pipelines and should not be used
for untrusted input.

## Competitor Landscape

The fair comparison set has multiple categories:

| Category | Tools | What they prove | What they do not prove |
| --- | --- | --- | --- |
| Command-line FASTQ processing | `fastp`, `seqkit`, `seqtk`, `samtools import` | End-user throughput and practical ecosystem relevance | Library embedding cost or zero-copy batch ergonomics |
| Rust parsing libraries | `seq_io`, `noodles-fastq`, `rust-bio` | Rust API ergonomics, parser behavior, safety tradeoffs | Full command-line workflow throughput |
| Compression libraries | `flate2`, `libdeflater`, `noodles-bgzf`, htslib/bgzip | Backend decompression/compression behavior | FASTQ parser and packed side-channel behavior |

Known published or maintained references:

- `fastp`: C++ FASTQ preprocessor with quality control, trimming, filtering,
  reporting, and multi-threaded design. It is a workflow competitor, not a
  parser-only peer. Reference: https://doi.org/10.1093/bioinformatics/bty560
- `seqtk`: Heng Li's lightweight C FASTA/Q toolkit. Reference:
  https://github.com/lh3/seqtk
- `seq_io`: Rust FASTA/FASTQ parser focused on high-performance borrowed
  records. It supports single-line FASTQ records and documents allocation-avoidant
  readers. Reference: https://docs.rs/seq_io/latest/seq_io/
- `noodles-bgzf`: Rust BGZF implementation from the noodles ecosystem.
  Reference: https://docs.rs/noodles-bgzf/latest/noodles_bgzf/

Current local Rust-library peer evidence is captured in:

- `docs/benchmarks/rust-peers/summary.md`: synthetic raw FASTQ parser-library
  comparison.
- `docs/benchmarks/rust-peers-drosophila-r1/summary.md`: Drosophila R1 raw
  FASTQ parser-library comparison.

Those rows currently show stable-default Microraptor close to `seq_io`, ahead of
`noodles-fastq` and `bio` on the local Drosophila R1 pass, and slightly behind
`seq_io`/`noodles-fastq` on the synthetic raw FASTQ fixture. The diagnostic peer
mode showed validation and accessor overhead are not the main difference; the
loss was dominated by newline discovery. The stronger claim should still be the
unified raw/gzip/BGZF streaming and packed side-channel framework, not universal
parser-only dominance over every Rust FASTQ reader.

Current local command-line/compression evidence is captured in:

- `docs/benchmarks/drosophila-1m/summary.md`: raw paired Drosophila Illumina
  rows with `seqkit`, `seqtk`, `samtools`, and `fastp` comparator timings.
- `docs/benchmarks/drosophila-compressed/summary.md`: Drosophila-derived gzip
  paired rows and an above-threshold combined BGZF input.
- `docs/benchmarks/drosophila-read-types/summary.md`: real Drosophila Illumina
  PE, PacBio CLR, and ONT FASTQ rows with command-line comparator timings.
- `docs/benchmarks/CORPUS.md`: the local biological corpus ladder, including
  larger Illumina paired-end inputs and ONT/PacBio read-type coverage from
  `~/Projects/Benchmarks`.

These are useful local artifacts, not universal claims. They should be treated
as evidence that the benchmark machinery can exercise real raw/gzip/BGZF
biological inputs and external tools on this machine.

## Ruthless Claim Boundaries

Do claim:

- Microraptor provides a compact Rust library core for raw/gzip/BGZF FASTQ
  streaming.
- It has a single parser path after decompression.
- It exposes borrowed records and optional packed base/quality side channels.
- It has benchmark scripts that check checksum parity across parser and trusted
  pack paths.
- Its benchmark suite includes both command-line workflow comparators and Rust
  parser-library peer comparators.

Do not claim without new evidence:

- Faster than `fastp`, `seqkit`, `seqtk`, `samtools`, `seq_io`, or noodles on
  all workloads.
- Faster on real biological datasets based only on synthetic fixtures.
- Suitable for multiline FASTQ.
- Suitable for reordered paired-end synchronization.
- Production-ready for every sequencing platform.

Required evidence before any public performance claim:

- Hardware, OS, compiler, feature flags, and command lines.
- Synthetic fixture results plus at least one real short-read R1/R2 gzip pair.
- At least one BGZF input large enough to cross the adaptive parallel threshold.
- Comparator versions and exact comparator commands.
- Checksums or equivalent record/base-count parity.
- A figure generated from committed scripts, not hand-built spreadsheets.
