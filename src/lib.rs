#![cfg_attr(feature = "simd", feature(portable_simd))]

pub mod benchutil;
#[cfg(feature = "bgzf")]
mod bgzf;
mod error;
mod fastq;
pub mod pack;
mod scan;
mod source;

#[cfg(feature = "bgzf")]
pub use bgzf::{
    BGZF_EOF_BLOCK, BgzfIndex, BgzfIndexEntry, BgzfInflateBackend, BgzfParallelReader, BgzfReader,
    BgzfVirtualOffset, BgzfWriter, build_bgzf_index, compress_bgzf_parallel,
    decompress_bgzf_parallel, decompress_bgzf_parallel_with_inflate_backend,
};
pub use error::{FastqError, FastqPosition, Result};
pub use fastq::{
    FastqBatch, FastqConfig, FastqPair, FastqReader, FastqRecord, InterleavedPairs, PairedRecords,
    PairingMode, RecordRef, paired_records, strip_pair_suffix,
};
pub use source::{open_fastq, open_fastq_with_config};
#[cfg(feature = "bgzf")]
pub use source::{open_fastq_bgzf_parallel, open_fastq_bgzf_parallel_with_config};
