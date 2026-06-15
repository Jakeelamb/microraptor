//! Optional memory-mapped resident visitors.

use std::fs::File;
use std::path::Path;

use memmap2::Mmap;

use crate::{
    FastaStats, FastaVisitRecord, FastqConfig, FastqVisitRecord, Result, count_fasta_bytes,
    visit_fasta_bytes, visit_fastq_bytes,
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
    visit_fasta_bytes(&map, visit)
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

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;

    use super::*;
    use crate::{FastqConfig, count_fasta_bytes, visit_fasta_bytes, visit_fastq_bytes};

    fn temp_file(name: &str, bytes: &[u8]) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "microraptor-mmap-{name}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::write(&path, bytes).unwrap();
        path
    }

    #[test]
    fn fastq_mmap_matches_resident_visitor() {
        let input = b"@r1\nACGT\n+\nIIII\n@r2\nTG\n+\n!!\n";
        let path = temp_file("fastq", input);
        let mut mmap_records = Vec::new();
        let mut resident_records = Vec::new();

        let mmap_count = visit_fastq_mmap(&path, FastqConfig::default(), |record| {
            mmap_records.push((record.name().to_vec(), record.seq().to_vec()));
            Ok(())
        })
        .unwrap();
        let resident_count = visit_fastq_bytes(input, FastqConfig::default(), |record| {
            resident_records.push((record.name().to_vec(), record.seq().to_vec()));
            Ok(())
        })
        .unwrap();

        assert_eq!(mmap_count, resident_count);
        assert_eq!(mmap_records, resident_records);
    }

    #[test]
    fn fasta_mmap_matches_resident_visitor_and_counter() {
        let input = b">seq1\nAC\nGT\n>seq2\nTTA\n";
        let path = temp_file("fasta", input);
        let mut mmap_records = Vec::new();
        let mut resident_records = Vec::new();

        let mmap_count = visit_fasta_mmap(&path, |record| {
            mmap_records.push((record.id_token().to_vec(), record.seq().to_vec()));
            Ok(())
        })
        .unwrap();
        let resident_count = visit_fasta_bytes(input, |record| {
            resident_records.push((record.id_token().to_vec(), record.seq().to_vec()));
            Ok(())
        })
        .unwrap();

        assert_eq!(mmap_count, resident_count);
        assert_eq!(mmap_records, resident_records);
        assert_eq!(
            count_fasta_mmap(&path).unwrap(),
            count_fasta_bytes(input).unwrap()
        );
    }
}
