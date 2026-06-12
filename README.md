# microraptor

Microraptor is a slab-based FASTQ streaming core. It is built around one invariant:
raw FASTQ, gzip FASTQ, and BGZF/ISA-L paths all produce the same
decompressed byte slabs, and the parser only sees bytes.

Current slice:

- raw FASTQ from any `Read`
- gzip FASTQ from paths with gzip magic, via `flate2::read::MultiGzDecoder`
- BGZF FASTQ detected by BGZF headers before ordinary gzip
- serial BGZF streaming reader and writer
- parallel BGZF decompression/compression entry points for independent blocks
- BGZF block index construction and `BgzfSeekReader` virtual-offset reads
- optional libdeflate BGZF inflate and deflate backends
- explicit buffered libdeflate gzip opener for bounded single gzip inputs
- reusable slab buffer with carry handling for records crossing slab boundaries
- SIMD newline scan on nightly through `std::simd`, with scalar fallback when the
  `simd` feature is disabled
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
- structured FASTQ parse errors with byte offset, record index, and line index
- zero-copy FASTQ record-id helpers for raw names, first tokens, and pair-normalized IDs
- stateful separate-file paired reader and interleaved FASTQ iterators with
  normalized-id validation
- configurable `PairValidation` modes for full ID checks, fast `/1` `/2`
  ordered-mate checks, or trusted ordered inputs
- minimal `FastqBatchSource` and `FastqPairBatchSource` traits for downstream modules

Backend boundary:

- ordinary gzip auto-open uses streaming `flate2`; `open_fastq_gzip_libdeflate`
  is explicit because it buffers the decompressed input
- BGZF is already block-aware and has parallel whole-input compression and
  decompression helpers, a bounded streaming parallel reader, an adaptive
  serial/parallel reader, and virtual-offset indexing/seek reads
- output compression supports flate2 by default and libdeflate when requested

Default streamer boundary:

- `open_fastq_with_config` is the default streamer surface for raw, gzip, and
  BGZF paths. Treat it as frozen for the current tiny-module scope unless a real
  workload exposes a correctness issue or measured bottleneck.
- performance work should preserve the one-streamer shape: scripts and benchmark
  rows may compare explicit alternatives, but production callers should not need
  to choose between competing default FASTQ pipelines.
- new streamer behavior needs parity evidence across default features,
  `--all-features`, and `--no-default-features`; timing-only wins are not enough
  to justify API churn.

Features:

- `simd`: nightly portable-SIMD newline scanner
- `gzip`: ordinary gzip input by gzip magic
- `bgzf`: BGZF reader, writer, detection, and parallel block helpers
- `libdeflate`: optional libdeflate BGZF inflate/deflate backends and explicit
  buffered gzip opener; makes BGZF auto-open use the fastest available inflate backend

Current limitations:

- FASTQ records must be four physical lines. Multiline sequence or quality
  fields are rejected as malformed or truncated input.
- Paired R1/R2 support validates ordered mates; it does not synchronize files
  with reordered records.

Benchmarking:

- `cargo bench --all-features`
- `cargo run --release --bin microraptor-bench -- --records 500000 --iters 7`
- `scripts/bench.sh`
- `scripts/benchmark-gauntlet.sh`
- `scripts/profile-perf.sh`
- `scripts/profile-hotpath.sh`

See `BENCHMARKING.md` for profiling details and result interpretation.

Robustness:

- `cargo fuzz run fastq_reader`
- `cargo fuzz run pack`
- `cargo fuzz run bgzf_roundtrip`

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
