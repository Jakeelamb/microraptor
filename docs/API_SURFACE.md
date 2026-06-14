# Public API Surface

This document records the intended public surface for `microraptor` `0.1.x`.
It is a release audit artifact, not a promise that the crate is frozen before
`1.0`. The purpose is to keep the first public release from accidentally
exporting implementation details that are hard to retract.

## Tier 1: Primary User Surface

These APIs are the default entry points for downstream scientific tools.

| API | Role | 0.1.x decision |
| --- | --- | --- |
| `FastqReader` | Borrowed batch reader over any `Read` source | Keep public |
| `FastqReader::visit_records`, `FastqVisitRecord` | Single-pass streaming visitor without batch side-table construction | Keep public |
| `visit_fastq_bytes` | Zero-copy visitor for complete resident FASTQ byte buffers | Keep public |
| `FastqBatch`, `FastqRecord`, `RecordRef` | Zero-copy record access within a reusable slab | Keep public |
| `FastqConfig` | Slab size, validation, and pairing configuration | Keep public |
| `open_fastq`, `open_fastq_with_config` | File-path opener for raw, gzip, and BGZF inputs | Keep public |
| `FastaReader`, `FastaBatch`, `FastaRecord`, `FastaRecordRef` | Streaming multiline FASTA batches over any `Read` source | Keep public |
| `FastaConfig` | FASTA batch, input-buffer, and sequence-length hint configuration | Keep public |
| `visit_fasta_bytes` | Zero-copy resident FASTA visitor with multiline folding fallback | Keep public |
| `visit_fasta_bytes_auto`, `detect_fasta_shape`, `FastaShape` | Resident FASTA shape detection and automatic strict two-line dispatch | Keep public |
| `visit_two_line_fasta_bytes`, `visit_two_line_fasta_read` | Strict `>header`/`sequence` fast paths for canonical two-line FASTA | Keep public |
| `count_two_line_fasta_bytes`, `count_two_line_fasta_read`, `FastaStats` | Strict two-line FASTA count/total-bases/light-checksum paths | Keep public |
| `FastaRecordSink`, `FastaVisitRecord` | Borrowed FASTA visitor sink and record view | Keep public |
| `open_fasta`, `open_fasta_with_config` | File-path opener for raw, gzip, and BGZF FASTA inputs | Keep public |
| `PairedFastqReader`, `PairedFastqBatch`, `FastqPair` | Ordered paired-end streaming | Keep public |
| `open_paired_fastq*` | Paired file opener variants | Keep public |
| `PairingMode`, `PairValidation`, `strip_pair_suffix` | Explicit ordered-pair validation semantics | Keep public |
| `FastqError`, `FastqPosition`, `Result` | Typed error reporting | Keep public |
| `FastqBatchSource`, `FastqPairBatchSource` | Generic adapters for downstream pipelines | Keep public |

Rationale: these are the crate's core value proposition. They expose borrowed
FASTQ and FASTA batches without imposing a workflow, allocator, or owned record
model.

Compression backends are not owned by microraptor. `flate2` and `libdeflate`
are third-party transport engines. Microraptor owns the parser APIs, BGZF
orchestration, backend-selection surface, and benchmark labels around those
engines.

## Tier 2: Advanced But Intentional Surface

These APIs are lower level, but they are still part of the intended scientific
systems surface because they let downstream tools avoid reparsing or
reallocating.

| API | Role | 0.1.x decision |
| --- | --- | --- |
| `pack::pack_bases`, `pack_bases_into`, `pack_bases_into_slices` | Two-bit base packing | Keep public |
| `pack::packed_base_at`, `is_masked` | Decode/access packed bases and ambiguity mask | Keep public |
| `pack::summarize_qualities`, `bin_qualities_into*` | Phred+33 summaries and bins | Keep public |
| `pack::PackedSequence`, `BaseSummary`, `QualitySummary`, `PackedRecordSummary` | Data carriers for packed side channels | Keep public |
| `pack::TrustedPackSink`, `TrustedPackedRecord`, `TrustedPackedPair`, `TrustedPackSlab` | Callback surface for streaming pack paths | Keep public with "trusted FASTQ" naming |
| `pack::pack_trusted_fastq*` | High-throughput four-line pack paths | Keep public with explicit trusted-input contract |
| `pack::selected_pack_kernel`, `PackKernel` | Build/host pack-kernel introspection | Keep public |

Rationale: the pack module is not a private optimization. It is a separable
side-channel API for tools that need compact base representation, ambiguity
masks, and quality summaries. The unsafe-looking part is semantic, not memory
unsafe: "trusted" means ordinary four-line FASTQ shape has already been chosen
as a workload contract. Public docs and benchmark text must keep that boundary
visible.

## Tier 3: Transport And Format Surface

These APIs make BGZF a first-class transport instead of hiding it behind file
auto-detection.

| API | Role | 0.1.x decision |
| --- | --- | --- |
| `BgzfReader`, `BgzfAutoReader`, `BgzfParallelReader` | Serial, adaptive, and bounded parallel BGZF decode | Keep public |
| `BgzfWriter`, `BGZF_EOF_BLOCK` | BGZF output support | Keep public |
| `BgzfInflateBackend`, `BgzfDeflateBackend` | Backend selection when `libdeflate` is enabled | Keep public |
| `BgzfParallelConfig`, `BgzfPipelineMetrics*` | Parallel threshold, queue, backend, and backpressure tuning | Keep public |
| `BgzfVirtualOffset`, `BgzfIndex`, `BgzfIndexEntry`, `build_bgzf_index` | Seek/index support | Keep public |
| `open_fastq_bgzf_*` | Explicit BGZF openers for benchmarking and tuning | Keep public |
| `compress_bgzf_parallel*`, `decompress_bgzf_parallel*` | Whole-buffer helpers for fixtures and controlled conversions | Keep public, but not the main streaming path |

Rationale: BGZF is common enough in bioinformatics that downstream users need
explicit knobs for backend choice, seek/index behavior, and bounded parallelism.
The public docs should continue to steer ordinary FASTQ consumers toward
`open_fastq` and advanced users toward explicit BGZF APIs only when they need
transport control.

## Hidden Or Internal Surface

`benchutil` is `#[doc(hidden)]` and exists for repository benches and generated
peer benchmark harnesses. It is not part of the documented user contract.

Internal parser scanners, slab framing helpers, and BGZF worker plumbing remain
private. New public exports should be added only when a downstream caller can
state a concrete use case that cannot be served by the existing batch, pack, or
transport tiers.

## Risks To Revisit Before 1.0

- `RecordRef` exposes raw byte ranges. This is valuable for zero-copy callers,
  but a future 1.0 API may prefer an opaque view if range layout changes.
- Trusted pack functions assume ordinary four-line FASTQ. They are useful and
  benchmarked, but public examples must not imply multiline FASTQ support.
- FASTA support is intentionally parser-only in `0.1.x`: no quality summaries,
  paired validation, or trusted FASTQ pack APIs apply to FASTA records.
- The `visit_two_line_fasta_*` and `count_two_line_fasta_*` functions are
  deliberately strict fast paths for canonical two-line FASTA. Use
  `FastaReader` or `visit_fasta_bytes` for ordinary multiline FASTA.
- Whole-buffer BGZF helpers are convenient for fixtures and conversions, but
  large production workflows should prefer streaming readers/writers.
- `PairValidation::FastSlash` is intentionally fast and narrow. Broader naming
  conventions should be modeled as new explicit validation modes, not hidden
  behavior changes.

## Release Decision

For `0.1.0`, keep the current public surface. It is broad, but the breadth maps
to real framework axes: parsing, ordered pairing, packed side channels, and
BGZF transport. The release should avoid any claim that all public APIs are
equally high level. README, rustdoc, and benchmark prose should keep directing
typical users to the Tier 1 APIs.
