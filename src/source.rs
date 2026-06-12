use std::fs::File;
use std::io::Read;
#[cfg(any(feature = "bgzf", feature = "gzip"))]
use std::io::{Seek, SeekFrom};
use std::path::Path;

use crate::error::Result;
use crate::fastq::{FastqConfig, FastqReader};
#[cfg(feature = "bgzf")]
use crate::{BgzfInflateBackend, BgzfParallelReader, BgzfReader, bgzf::is_bgzf_header};

pub fn open_fastq(path: impl AsRef<Path>) -> Result<FastqReader<Box<dyn Read + Send>>> {
    open_fastq_with_config(path, FastqConfig::default())
}

pub fn open_fastq_with_config(
    path: impl AsRef<Path>,
    config: FastqConfig,
) -> Result<FastqReader<Box<dyn Read + Send>>> {
    #[cfg(any(feature = "bgzf", feature = "gzip"))]
    let mut file = File::open(path)?;
    #[cfg(not(any(feature = "bgzf", feature = "gzip")))]
    let file = File::open(path)?;

    #[cfg(any(feature = "bgzf", feature = "gzip"))]
    {
        let mut prefix = [0_u8; 18];
        let n = file.read(&mut prefix)?;
        file.seek(SeekFrom::Start(0))?;

        #[cfg(feature = "bgzf")]
        if is_bgzf_header(&prefix[..n]) {
            let reader: Box<dyn Read + Send> = Box::new(BgzfReader::new(file));
            return Ok(FastqReader::with_config(reader, config));
        }

        #[cfg(feature = "gzip")]
        if n >= 2 && prefix[..2] == [0x1f, 0x8b] {
            let reader: Box<dyn Read + Send> = Box::new(flate2::read::MultiGzDecoder::new(file));
            return Ok(FastqReader::with_config(reader, config));
        }
    }

    let reader: Box<dyn Read + Send> = Box::new(file);
    Ok(FastqReader::with_config(reader, config))
}

#[cfg(feature = "bgzf")]
pub fn open_fastq_bgzf_parallel(
    path: impl AsRef<Path>,
    workers: usize,
) -> Result<FastqReader<BgzfParallelReader>> {
    open_fastq_bgzf_parallel_with_config(path, workers, FastqConfig::default())
}

#[cfg(feature = "bgzf")]
pub fn open_fastq_bgzf_parallel_with_config(
    path: impl AsRef<Path>,
    workers: usize,
    config: FastqConfig,
) -> Result<FastqReader<BgzfParallelReader>> {
    open_fastq_bgzf_parallel_with_backend(path, workers, BgzfInflateBackend::default(), config)
}

#[cfg(feature = "bgzf")]
pub fn open_fastq_bgzf_flate2(path: impl AsRef<Path>) -> Result<FastqReader<BgzfReader<File>>> {
    open_fastq_bgzf_with_backend(path, BgzfInflateBackend::Flate2, FastqConfig::default())
}

#[cfg(feature = "bgzf")]
pub fn open_fastq_bgzf_with_backend(
    path: impl AsRef<Path>,
    backend: BgzfInflateBackend,
    config: FastqConfig,
) -> Result<FastqReader<BgzfReader<File>>> {
    let file = File::open(path)?;
    Ok(FastqReader::with_config(
        BgzfReader::with_inflate_backend(file, backend),
        config,
    ))
}

#[cfg(feature = "bgzf")]
pub fn open_fastq_bgzf_parallel_with_backend(
    path: impl AsRef<Path>,
    workers: usize,
    backend: BgzfInflateBackend,
    config: FastqConfig,
) -> Result<FastqReader<BgzfParallelReader>> {
    let file = File::open(path)?;
    Ok(FastqReader::with_config(
        BgzfParallelReader::with_inflate_backend(file, workers, backend)?,
        config,
    ))
}

#[cfg(all(test, feature = "gzip"))]
mod tests {
    use std::io::Write;

    use super::*;

    #[test]
    #[cfg(feature = "gzip")]
    fn opens_gzip_by_magic() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("microraptor-{}.fq.gz", std::process::id()));
        let file = File::create(&path).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::fast());
        encoder.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
        encoder.finish().unwrap();

        let mut reader = open_fastq(&path).unwrap();
        let batch = reader.next_batch().unwrap().unwrap();
        let rec = batch.records().next().unwrap();
        assert_eq!(rec.seq(), b"ACGT");

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    #[cfg(feature = "bgzf")]
    fn opens_bgzf_by_magic() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("microraptor-{}.bgz", std::process::id()));
        let mut writer = crate::BgzfWriter::new(File::create(&path).unwrap());
        writer.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
        writer.finish().unwrap();

        let mut reader = open_fastq(&path).unwrap();
        let batch = reader.next_batch().unwrap().unwrap();
        let rec = batch.records().next().unwrap();
        assert_eq!(rec.seq(), b"ACGT");

        let mut parallel = open_fastq_bgzf_parallel(&path, 2).unwrap();
        let batch = parallel.next_batch().unwrap().unwrap();
        let rec = batch.records().next().unwrap();
        assert_eq!(rec.seq(), b"ACGT");

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    #[cfg(all(feature = "bgzf", feature = "libdeflate"))]
    fn explicit_bgzf_backend_openers_parse_same_file() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("microraptor-backend-{}.bgz", std::process::id()));
        let mut writer = crate::BgzfWriter::new(File::create(&path).unwrap());
        writer
            .write_all(b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nJJJJ\n")
            .unwrap();
        writer.finish().unwrap();

        let mut auto = open_fastq(&path).unwrap();
        let mut flate2 = open_fastq_bgzf_flate2(&path).unwrap();
        let mut libdeflate = open_fastq_bgzf_with_backend(
            &path,
            BgzfInflateBackend::Libdeflate,
            FastqConfig::default(),
        )
        .unwrap();

        let auto_stats = crate::benchutil::consume_fastq(&mut auto).unwrap();
        let flate2_stats = crate::benchutil::consume_fastq(&mut flate2).unwrap();
        let libdeflate_stats = crate::benchutil::consume_fastq(&mut libdeflate).unwrap();

        assert_eq!(auto_stats, flate2_stats);
        assert_eq!(auto_stats, libdeflate_stats);

        std::fs::remove_file(path).unwrap();
    }
}
