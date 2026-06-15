use std::io::Read;

use crate::error::Result;
use crate::fasta::{FastaBatch, FastaReader};
use crate::fastq::{FastqBatch, FastqReader, PairedFastqBatch, PairedFastqReader};

/// Minimal abstraction over a source of FASTA batches.
///
/// Downstream pipeline stages can depend on this trait instead of the concrete
/// [`FastaReader`] type when they only need batch iteration.
pub trait FastaBatchSource {
    /// Read the next FASTA batch.
    fn next_fasta_batch(&mut self) -> Result<Option<FastaBatch<'_>>>;
}

impl<R: Read> FastaBatchSource for FastaReader<R> {
    fn next_fasta_batch(&mut self) -> Result<Option<FastaBatch<'_>>> {
        self.next_batch()
    }
}

/// Minimal abstraction over a source of single-end FASTQ batches.
///
/// Downstream pipeline stages can depend on this trait instead of the concrete
/// [`FastqReader`] type when they only need batch iteration.
pub trait FastqBatchSource {
    /// Read the next single-end batch.
    fn next_fastq_batch(&mut self) -> Result<Option<FastqBatch<'_>>>;
}

impl<R: Read> FastqBatchSource for FastqReader<R> {
    fn next_fastq_batch(&mut self) -> Result<Option<FastqBatch<'_>>> {
        self.next_batch()
    }
}

/// Minimal abstraction over a source of ordered paired-end FASTQ batches.
pub trait FastqPairBatchSource {
    /// Read the next paired-end batch.
    fn next_fastq_pair_batch(&mut self) -> Result<Option<PairedFastqBatch<'_>>>;
}

impl<R1: Read, R2: Read> FastqPairBatchSource for PairedFastqReader<R1, R2> {
    fn next_fastq_pair_batch(&mut self) -> Result<Option<PairedFastqBatch<'_>>> {
        self.next_pair_batch()
    }
}
