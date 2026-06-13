#![cfg_attr(feature = "simd", feature(portable_simd))]
#![warn(missing_docs)]
//! Slab-based FASTQ streaming and packing for raw, gzip, and BGZF inputs.
//!
//! Microraptor is a library core for downstream scientific tools that need a
//! low-allocation stream of ordinary short-read FASTQ records. Decompression is
//! treated as a transport layer: raw FASTQ, gzip FASTQ, and BGZF FASTQ are
//! converted into byte slabs, and the parser consumes the same slab model for
//! all input backends.
//!
//! The default feature set builds on stable Rust and includes gzip and BGZF
//! input support. Nightly-only SIMD acceleration is available through the
//! explicit `simd` feature.
//!
//! # Choosing an entry point
//!
//! Use [`FastqReader`] when you already have an object that implements
//! [`std::io::Read`]. Use [`visit_fastq_bytes`] when a complete FASTQ byte
//! buffer is already resident in memory. Use [`open_fastq`] or
//! [`open_fastq_with_config`] when you want file-path auto-detection for raw
//! FASTQ, ordinary gzip, and BGZF. Use [`PairedFastqReader`] or
//! [`open_paired_fastq`] for ordered R1/R2 streams.
//!
//! # Scope
//!
//! Microraptor parses four-line FASTQ records. It validates ordered paired-end
//! reads in separate R1/R2 streams or adjacent interleaved records, but it does
//! not synchronize reordered mates. It does not trim adapters, filter reads,
//! align reads, or generate quality-control reports.
//!
//! # Lifetimes and allocation
//!
//! Batches borrow from the reader's reusable slab buffer. A [`FastqBatch`] is
//! valid until the next mutable reader call. Clone or copy record data if it
//! must outlive the batch. This design keeps the parser low-allocation, but it
//! means callers should process each batch before advancing the reader.
//!
//! # Feature flags
//!
//! - `gzip` enables ordinary gzip auto-detection and streaming decode.
//! - `bgzf` enables BGZF readers, writers, indexing, and adaptive parallel
//!   decoding.
//! - `libdeflate` enables optional libdeflate BGZF backends and an explicit
//!   buffered gzip opener.
//! - `simd` enables nightly portable-SIMD scanner and packing paths.
//!
//! # Example
//!
//! ```
//! use microraptor::FastqReader;
//!
//! let data = b"@r1\nACGT\n+\nIIII\n";
//! let mut reader = FastqReader::new(&data[..]);
//! let mut records = 0;
//!
//! while let Some(batch) = reader.next_batch()? {
//!     for record in batch.records() {
//!         assert_eq!(record.seq(), b"ACGT");
//!         records += 1;
//!     }
//! }
//!
//! assert_eq!(records, 1);
//! # Ok::<(), microraptor::FastqError>(())
//! ```
//!
//! # Paired reads
//!
//! ```
//! use microraptor::PairedFastqReader;
//!
//! let r1 = b"@frag/1\nACGT\n+\nIIII\n";
//! let r2 = b"@frag/2\nTGCA\n+\nJJJJ\n";
//! let mut reader = PairedFastqReader::new(&r1[..], &r2[..]);
//! let batch = reader.next_pair_batch()?.expect("one paired batch");
//! let pair = batch.pairs().next().expect("one read pair");
//!
//! assert_eq!(pair.pair_id(), b"frag");
//! assert_eq!(pair.first().seq(), b"ACGT");
//! assert_eq!(pair.second().seq(), b"TGCA");
//! # Ok::<(), microraptor::FastqError>(())
//! ```

#[doc(hidden)]
pub mod benchutil;
#[cfg(feature = "bgzf")]
mod bgzf;
mod error;
mod fastq;
/// Base/quality packing and trusted four-line FASTQ pack paths.
///
/// The high-level FASTQ readers expose borrowed records. This module provides
/// the lower-level side-channel representation used when downstream code wants
/// packed two-bit bases, ambiguity masks, and Phred+33 summaries without
/// allocating an owned record per read.
pub mod pack;
mod scan;
mod source;
mod stream;

#[cfg(feature = "bgzf")]
pub use bgzf::{
    BGZF_EOF_BLOCK, BgzfAutoReader, BgzfDeflateBackend, BgzfIndex, BgzfIndexEntry,
    BgzfInflateBackend, BgzfParallelConfig, BgzfParallelReader, BgzfPipelineMetrics,
    BgzfPipelineMetricsSnapshot, BgzfReader, BgzfSeekReader, BgzfVirtualOffset, BgzfWriter,
    build_bgzf_index, compress_bgzf_parallel, compress_bgzf_parallel_with_deflate_backend,
    decompress_bgzf_parallel, decompress_bgzf_parallel_with_inflate_backend,
};
pub use error::{FastqError, FastqPosition, Result};
pub use fastq::{
    FastqBatch, FastqConfig, FastqPair, FastqReader, FastqRecord, FastqVisitRecord,
    InterleavedPairs, PairValidation, PairedFastqBatch, PairedFastqPairs, PairedFastqReader,
    PairedRecords, PairingMode, RecordRef, paired_records, strip_pair_suffix, visit_fastq_bytes,
};
pub use source::{
    open_fastq, open_fastq_with_config, open_paired_fastq, open_paired_fastq_with_config,
    open_paired_fastq_with_configs,
};
#[cfg(feature = "bgzf")]
pub use source::{
    open_fastq_bgzf_adaptive, open_fastq_bgzf_flate2, open_fastq_bgzf_parallel,
    open_fastq_bgzf_parallel_with_backend, open_fastq_bgzf_parallel_with_config,
    open_fastq_bgzf_parallel_with_options, open_fastq_bgzf_with_backend,
};
#[cfg(all(feature = "gzip", feature = "libdeflate"))]
pub use source::{open_fastq_gzip_libdeflate, open_fastq_gzip_libdeflate_with_config};
pub use stream::{FastqBatchSource, FastqPairBatchSource};
