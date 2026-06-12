#![cfg_attr(feature = "simd", feature(portable_simd))]

pub mod benchutil;
#[cfg(feature = "bgzf")]
mod bgzf;
mod error;
mod fastq;
pub mod pack;
mod scan;
mod source;
mod stream;

#[cfg(feature = "bgzf")]
pub use bgzf::{
    BGZF_EOF_BLOCK, BgzfDeflateBackend, BgzfIndex, BgzfIndexEntry, BgzfInflateBackend,
    BgzfParallelConfig, BgzfParallelReader, BgzfReader, BgzfSeekReader, BgzfVirtualOffset,
    BgzfWriter, build_bgzf_index, compress_bgzf_parallel,
    compress_bgzf_parallel_with_deflate_backend, decompress_bgzf_parallel,
    decompress_bgzf_parallel_with_inflate_backend,
};
pub use error::{FastqError, FastqPosition, Result};
pub use fastq::{
    FastqBatch, FastqConfig, FastqPair, FastqReader, FastqRecord, InterleavedPairs,
    PairedFastqBatch, PairedFastqPairs, PairedFastqReader, PairedRecords, PairingMode, RecordRef,
    paired_records, strip_pair_suffix,
};
pub use source::{
    open_fastq, open_fastq_with_config, open_paired_fastq, open_paired_fastq_with_config,
    open_paired_fastq_with_configs,
};
#[cfg(feature = "bgzf")]
pub use source::{
    open_fastq_bgzf_flate2, open_fastq_bgzf_parallel, open_fastq_bgzf_parallel_with_backend,
    open_fastq_bgzf_parallel_with_config, open_fastq_bgzf_parallel_with_options,
    open_fastq_bgzf_with_backend,
};
#[cfg(all(feature = "gzip", feature = "libdeflate"))]
pub use source::{open_fastq_gzip_libdeflate, open_fastq_gzip_libdeflate_with_config};
pub use stream::{FastqBatchSource, FastqPairBatchSource};
