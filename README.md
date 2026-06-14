# microraptor

Microraptor is a streaming FASTQ and FASTA parsing core for Rust. It is designed
for downstream scientific tools that need raw, gzip, and BGZF sequence inputs to
converge on low-allocation parser paths.

The framework is built around one invariant: decompression is a transport layer.
Raw, gzip, and BGZF paths produce the same decompressed byte stream for each
format-specific parser.

Compression engines are third-party backends, not microraptor inventions.
`flate2` provides the default Rust gzip/BGZF transport path. `libdeflate` is the
upstream high-performance DEFLATE implementation used through the Rust
`libdeflater` wrapper when the optional `libdeflate` feature or a benchmark row
explicitly selects it. Microraptor's own code is the FASTQ/FASTA parser surface,
BGZF orchestration, backend selection, and benchmark harness around those
compression libraries.

Microraptor is intentionally not an all-in-one preprocessing suite. It does not
trim adapters, filter reads, produce QC reports, align reads, or synchronize
reordered paired-end files. Its useful surface is narrower: validated FASTQ
batches, streaming multiline FASTA batches, ordered mate validation, BGZF-aware
input, and optional FASTQ packed base/quality side channels.

Publication surfaces:

- [docs/FRAMEWORK.md](docs/FRAMEWORK.md): framework analysis, competitors, and
  claim boundaries.
- [BENCHMARKING.md](BENCHMARKING.md): benchmark commands, interpretation, and
  reproducibility protocol.
- [docs/API_SURFACE.md](docs/API_SURFACE.md): public API tiering and `0.1.x`
  release-surface decision.
- [docs/REPLICATION.md](docs/REPLICATION.md): release-commit regeneration,
  replication-kit, and independent-machine evidence protocol.
- [docs/PUBLISHING.md](docs/PUBLISHING.md): release-readiness checklist.
- [docs/RELEASE_AUDIT.md](docs/RELEASE_AUDIT.md): current requirement matrix
  and remaining public-release gaps.
- [CHANGELOG.md](CHANGELOG.md): release notes and benchmark-evidence summary.
- [SECURITY.md](SECURITY.md): vulnerability reporting scope for parser, pack,
  and BGZF issues.
- [CITATION.cff](CITATION.cff): citation metadata for scientific use.

Current slice:

- raw FASTQ from any `Read`
- gzip FASTQ from paths with gzip magic, via `flate2::read::MultiGzDecoder`
- BGZF FASTQ detected by BGZF headers before ordinary gzip
- raw/gzip/BGZF FASTA via `FastaReader`, `open_fasta`, and
  `open_fasta_with_config`
- resident FASTA visitors, including strict two-line FASTA fast paths for
  canonical `>header`/`sequence` inputs
- light-accounting FASTA stats and strict resident/streaming two-line counters
  for count/total-bases/checksum workloads
- serial BGZF streaming reader and writer
- parallel BGZF decompression/compression entry points for independent blocks
- BGZF block index construction and `BgzfSeekReader` virtual-offset reads
- optional libdeflate BGZF inflate and deflate backends
- explicit buffered libdeflate gzip opener for bounded single gzip inputs
- reusable slab buffer with carry handling for records crossing slab boundaries
- stable default `memchr` newline scanner, with nightly `std::simd`
  acceleration behind the explicit `simd` feature
- borrowed `RecordRef` ranges instead of per-record allocation
- optional 2-bit base packing with an ambiguity mask
- Phred+33 quality summaries and threshold binning
- trusted streaming FASTQ pack path for four-line records without batch
  allocation
- trusted packed-record sink APIs, paired trusted packing, and selected pack
  kernel reporting
- direct single-pass trusted pack scanner and assembly audit script for pack
  kernel inspection
- lockstep streaming paired trusted pack path with bounded mate buffering
- fused base+quality packing with compact quad LUTs, AVX2 quality reductions,
  and slab/BGZF pack benchmark gates
- canonical A/C/G/T SIMD chunk packing fast path and concrete trusted stats sink
  for low-overhead pack benchmarks
- structured FASTQ/FASTA parse errors with byte offset, record index, and line index
- zero-copy FASTQ record-id helpers for raw names, first tokens, and pair-normalized IDs
- stateful separate-file paired reader and interleaved FASTQ iterators with
  normalized-id validation
- configurable `PairValidation` modes for full ID checks, fast `/1` `/2`
  ordered-mate checks, or trusted ordered inputs
- minimal `FastqBatchSource` and `FastqPairBatchSource` traits for downstream modules

Backend boundary:

- ordinary gzip auto-open uses streaming `flate2`, a third-party Rust
  compression crate; `open_fastq_gzip_libdeflate` is explicit because it buffers
  the decompressed input through the third-party `libdeflate` engine via the
  `libdeflater` Rust wrapper
- BGZF is already block-aware and has parallel whole-input compression and
  decompression helpers, a bounded streaming parallel reader, an adaptive
  serial/parallel reader, and virtual-offset indexing/seek reads
- output compression supports third-party `flate2` by default and third-party
  `libdeflate` when requested

Default streamer boundary:

- `open_fastq_with_config` is the default FASTQ streamer surface for raw, gzip,
  and BGZF paths. `open_fasta_with_config` is the matching FASTA opener for the
  same transport layer. Treat them as frozen for the current tiny-module scope
  unless a real workload exposes a correctness issue or measured bottleneck.
- performance work should preserve the one-streamer shape: scripts and benchmark
  rows may compare explicit alternatives, but production callers should not need
  to choose between competing default FASTQ pipelines.
- new streamer behavior needs parity evidence across default features,
  `--all-features`, and `--no-default-features`; timing-only wins are not enough
  to justify API churn.

Features:

- default: stable Rust raw/gzip/BGZF streaming
- `simd`: nightly portable-SIMD newline and pack-path acceleration
- `gzip`: ordinary gzip input by gzip magic
- `bgzf`: BGZF reader, writer, detection, and parallel block helpers
- `libdeflate`: optional libdeflate BGZF inflate/deflate backends and explicit
  buffered gzip opener through the `libdeflater` wrapper; makes BGZF auto-open
  use the fastest available configured inflate backend

Current limitations:

- FASTQ records must be four physical lines. Multiline sequence or quality
  fields are rejected as malformed or truncated input.
- FASTA sequence lines may be multiline and are concatenated per record. FASTA
  has no quality stream, so FASTQ quality summaries and trusted pack paths do
  not apply to FASTA records.
- `visit_two_line_fasta_bytes` and `visit_two_line_fasta_read` are strict fast
  paths for canonical two-line FASTA. Use `FastaReader` or `visit_fasta_bytes`
  when records may contain multiline sequences or blank lines.
- `count_two_line_fasta_bytes` and `count_two_line_fasta_read` are stricter
  count/stat paths for canonical two-line FASTA only. They do not replace the
  robust multiline parser.
- Paired R1/R2 support validates ordered mates; it does not synchronize files
  with reordered records.

Benchmarking:

- `cargo +nightly bench --all-features`
- `cargo run --release --bin microraptor-bench -- --records 500000 --iters 7`
- `scripts/bench.sh`
- `scripts/benchmark-gauntlet.sh`
- `scripts/render-benchmark-report.sh`
- `scripts/check-benchmark-snapshots.sh`
- `scripts/benchmark-rust-peers.sh`
- `scripts/benchmark-fasta-peers.sh`
- `scripts/check-replication-host.sh`
- `scripts/discover-local-benchmark-corpus.sh`
- `scripts/profile-perf.sh`
- `scripts/profile-hotpath.sh`
- `scripts/benchmark-rust-peer-size-sweep.sh`
- `scripts/benchmark-fasta-peer-size-sweep.sh`
- `scripts/benchmark-fasta-gauntlet.sh`
- `scripts/release-gate.sh`
- `scripts/export-replication-kit.sh`

See `BENCHMARKING.md` for profiling details and result interpretation.

Current local benchmark snapshots:

- [docs/benchmarks/drosophila-1m/summary.md](docs/benchmarks/drosophila-1m/summary.md):
  1M ordered Drosophila Illumina read pairs from the local benchmark corpus,
  plus installed `seqkit`, `seqtk`, `samtools`, and `fastp` wall-clock rows.
- [docs/benchmarks/drosophila-compressed/summary.md](docs/benchmarks/drosophila-compressed/summary.md):
  real Drosophila-derived gzip and BGZF rows, including a combined BGZF input
  above the adaptive parallel threshold.
- [docs/benchmarks/drosophila-read-types/summary.md](docs/benchmarks/drosophila-read-types/summary.md):
  real Drosophila Illumina PE, PacBio CLR, and ONT FASTQ rows from the local
  benchmark corpus with installed command-line comparator timings.
- [docs/benchmarks/fasta-peer-size-sweep/summary.md](docs/benchmarks/fasta-peer-size-sweep/summary.md):
  synthetic two-line FASTA raw/gzip parser-framework size sweep with Rust
  parser peers, an 8-thread cap, and a required microraptor-winner guard.
- [docs/benchmarks/fasta-gauntlet/summary.md](docs/benchmarks/fasta-gauntlet/summary.md):
  FASTA shape/transport gauntlet covering two-line DNA, wrapped DNA, many tiny
  records, long contigs, protein FASTA, raw/gzip/BGZF transport, memory smoke
  rows, and installed command-line comparator timings.
- [docs/benchmarks/independent-organisms/summary.md](docs/benchmarks/independent-organisms/summary.md):
  independent non-Drosophila E. coli and yeast paired FASTQ rows from the local
  benchmark corpus with installed command-line comparator timings.
- [docs/benchmarks/latest/summary.md](docs/benchmarks/latest/summary.md):
  deterministic synthetic fixture gauntlet.
- [docs/benchmarks/rust-peers/summary.md](docs/benchmarks/rust-peers/summary.md):
  synthetic raw FASTQ parser-library comparison against `seq_io`,
  `noodles-fastq`, and `bio`, including separate Microraptor streaming and
  resident-slice visitor rows.
- [docs/benchmarks/rust-peers-drosophila-r1/summary.md](docs/benchmarks/rust-peers-drosophila-r1/summary.md):
  Drosophila R1 raw FASTQ parser-library comparison against the same Rust
  peers.

For public performance claims, regenerate the gauntlet from the current commit,
render the benchmark summary/figure, and record hardware plus comparator
versions. Target artifacts checked into a local workspace are development
evidence, not publication evidence.

Robustness:

- `cargo fuzz run fastq_reader`
- `cargo fuzz run pack`
- `cargo fuzz run bgzf_roundtrip`

Citation:

If you use microraptor in scientific work, cite the repository and exact version
or commit used. Citation metadata is provided in [CITATION.cff](CITATION.cff).

Contributing:

See [CONTRIBUTING.md](CONTRIBUTING.md). Performance or benchmark contributions
must include exact commands, environment, comparator versions, raw outputs, and
generated figures; timing-only claims without reproducible evidence should not
be merged.

Example:

```rust
use microraptor::FastqReader;

let data = b"@r1\nACGT\n+\nIIII\n";
let mut reader = FastqReader::new(&data[..]);
while let Some(batch) = reader.next_batch()? {
    for record in batch.records() {
        assert_eq!(record.seq(), b"ACGT");
    }
}
# Ok::<(), microraptor::FastqError>(())
```

Separate R1/R2 streams can be read as stateful paired batches:

```rust
use microraptor::PairedFastqReader;

let r1 = b"@frag/1\nACGT\n+\nIIII\n";
let r2 = b"@frag/2\nTGCA\n+\nJJJJ\n";
let mut reader = PairedFastqReader::new(&r1[..], &r2[..]);
let batch = reader.next_pair_batch()?.unwrap();
for pair in batch.pairs() {
    assert_eq!(pair.pair_id(), b"frag");
}
# Ok::<(), microraptor::FastqError>(())
```

Interleaved paired-end batches can be enabled without changing the streaming
path:

```rust
use microraptor::{FastqConfig, FastqReader};

let data = b"@frag/1\nACGT\n+\nIIII\n@frag/2\nTGCA\n+\nJJJJ\n";
let mut reader = FastqReader::with_config(&data[..], FastqConfig::default().interleaved());
let batch = reader.next_batch()?.unwrap();
for pair in batch.interleaved_pairs()? {
    assert_eq!(pair.pair_id(), b"frag");
}
# Ok::<(), microraptor::FastqError>(())
```

FASTA streams use a separate reader:

```rust
use microraptor::FastaReader;

let data = b">seq1 description\nACG\nTN\n";
let mut reader = FastaReader::new(&data[..]);
let batch = reader.next_batch()?.unwrap();
for record in batch.records() {
    assert_eq!(record.id_token(), b"seq1");
    assert_eq!(record.seq(), b"ACGTN");
}
# Ok::<(), microraptor::FastqError>(())
```

## License

Licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.
