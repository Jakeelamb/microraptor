# microraptor

Microraptor is a slab-based FASTQ streaming core. It is built around one invariant:
raw FASTQ, gzip FASTQ, and future BGZF/ISA-L paths all produce the same
decompressed byte slabs, and the parser only sees bytes.

Current slice:

- raw FASTQ from any `Read`
- gzip FASTQ from paths with gzip magic, via `flate2::read::MultiGzDecoder`
- BGZF FASTQ detected by BGZF headers before ordinary gzip
- serial BGZF streaming reader and writer
- parallel BGZF decompression/compression entry points for independent blocks
- reusable slab buffer with carry handling for records crossing slab boundaries
- SIMD newline scan on nightly through `std::simd`, with scalar fallback when the
  `simd` feature is disabled
- borrowed `RecordRef` ranges instead of per-record allocation
- optional 2-bit base packing with an ambiguity mask
- Phred+33 quality summaries and threshold binning
- structured FASTQ parse errors with byte offset, record index, and line index
- zero-copy FASTQ record-id helpers for raw names, first tokens, and pair-normalized IDs
- paired and interleaved FASTQ iterators with normalized-id validation

Backend boundary:

- ordinary gzip currently uses `flate2`; the decoder stage is isolated so ISA-L
  can replace it without touching FASTQ framing
- BGZF is already block-aware and has parallel whole-input compression and
  decompression helpers, plus a bounded streaming parallel reader
- output compression currently uses Rust deflate through `flate2`; libdeflate or
  ISA-L can replace block compression later

Features:

- `simd`: nightly portable-SIMD newline scanner
- `gzip`: ordinary gzip input by gzip magic
- `bgzf`: BGZF reader, writer, detection, and parallel block helpers

Benchmarking:

- `cargo bench --all-features`
- `cargo run --release --bin microraptor-bench -- --records 500000 --iters 7`
- `scripts/profile-perf.sh`

See `BENCHMARKING.md` for profiling details and result interpretation.

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
