//! Optional memory-mapped resident visitors.

use std::fs::File;
use std::path::Path;

use memmap2::Mmap;

use crate::{
    FastaStats, FastaVisitRecord, FastqConfig, FastqVisitRecord, Result, count_fasta_bytes,
    visit_fasta_bytes_auto, visit_fastq_bytes,
};

/// Visit FASTQ records from a memory-mapped file.
///
/// This is intended for already-resident local files where mapping avoids an
/// extra read buffer and lets callers use the resident byte-slice parser.
pub fn visit_fastq_mmap<F>(path: impl AsRef<Path>, config: FastqConfig, visit: F) -> Result<u64>
where
    F: FnMut(FastqVisitRecord<'_>) -> Result<()>,
{
    let map = map_file(path)?;
    visit_fastq_bytes(&map, config, visit)
}

/// Visit FASTA records from a memory-mapped file.
pub fn visit_fasta_mmap<F>(path: impl AsRef<Path>, visit: F) -> Result<u64>
where
    F: FnMut(FastaVisitRecord<'_>) -> Result<()>,
{
    let map = map_file(path)?;
    visit_fasta_bytes_auto(&map, visit)
}

/// Count FASTA records and bases from a memory-mapped file.
pub fn count_fasta_mmap(path: impl AsRef<Path>) -> Result<FastaStats> {
    let map = map_file(path)?;
    count_fasta_bytes(&map)
}

fn map_file(path: impl AsRef<Path>) -> Result<Mmap> {
    let file = File::open(path)?;
    // SAFETY: The map is read-only and kept alive for the duration of parsing.
    // Callers must not mutate the file concurrently through another handle.
    unsafe { Mmap::map(&file) }.map_err(Into::into)
}
