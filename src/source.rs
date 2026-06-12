use std::fs::File;
#[cfg(all(feature = "gzip", feature = "libdeflate"))]
use std::io::Cursor;
use std::io::Read;
#[cfg(any(feature = "bgzf", feature = "gzip"))]
use std::io::{Seek, SeekFrom};
use std::path::Path;

use crate::error::Result;
use crate::fastq::{FastqConfig, FastqReader, PairedFastqReader};
#[cfg(feature = "bgzf")]
use crate::{
    BgzfAutoReader, BgzfInflateBackend, BgzfParallelConfig, BgzfParallelReader, BgzfReader,
    bgzf::is_bgzf_header,
};

#[cfg(any(feature = "bgzf", feature = "gzip"))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum InputKind {
    Raw,
    #[cfg(feature = "gzip")]
    Gzip,
    #[cfg(feature = "bgzf")]
    Bgzf,
}

pub fn open_fastq(path: impl AsRef<Path>) -> Result<FastqReader<Box<dyn Read + Send>>> {
    open_fastq_with_config(path, FastqConfig::default())
}

pub fn open_fastq_with_config(
    path: impl AsRef<Path>,
    config: FastqConfig,
) -> Result<FastqReader<Box<dyn Read + Send>>> {
    let reader = open_read_by_magic(path)?;
    Ok(FastqReader::with_config(reader, config))
}

pub fn open_paired_fastq(
    first_path: impl AsRef<Path>,
    second_path: impl AsRef<Path>,
) -> Result<PairedFastqReader<Box<dyn Read + Send>, Box<dyn Read + Send>>> {
    open_paired_fastq_with_config(first_path, second_path, FastqConfig::default())
}

pub fn open_paired_fastq_with_config(
    first_path: impl AsRef<Path>,
    second_path: impl AsRef<Path>,
    config: FastqConfig,
) -> Result<PairedFastqReader<Box<dyn Read + Send>, Box<dyn Read + Send>>> {
    open_paired_fastq_with_configs(first_path, config.clone(), second_path, config)
}

pub fn open_paired_fastq_with_configs(
    first_path: impl AsRef<Path>,
    first_config: FastqConfig,
    second_path: impl AsRef<Path>,
    second_config: FastqConfig,
) -> Result<PairedFastqReader<Box<dyn Read + Send>, Box<dyn Read + Send>>> {
    let first = open_read_by_magic(first_path)?;
    let second = open_read_by_magic(second_path)?;
    Ok(PairedFastqReader::with_configs(
        first,
        first_config,
        second,
        second_config,
    ))
}

fn open_read_by_magic(path: impl AsRef<Path>) -> Result<Box<dyn Read + Send>> {
    #[cfg(any(feature = "bgzf", feature = "gzip"))]
    let mut file = File::open(path)?;
    #[cfg(not(any(feature = "bgzf", feature = "gzip")))]
    let file = File::open(path)?;

    #[cfg(any(feature = "bgzf", feature = "gzip"))]
    {
        let kind = detect_input_kind(&mut file)?;
        file.seek(SeekFrom::Start(0))?;

        match kind {
            #[cfg(feature = "bgzf")]
            InputKind::Bgzf => return Ok(Box::new(BgzfReader::new(file))),
            #[cfg(feature = "gzip")]
            InputKind::Gzip => return Ok(Box::new(flate2::read::MultiGzDecoder::new(file))),
            InputKind::Raw => {}
        }
    }

    let reader: Box<dyn Read + Send> = Box::new(file);
    Ok(reader)
}

#[cfg(any(feature = "bgzf", feature = "gzip"))]
fn detect_input_kind(file: &mut File) -> Result<InputKind> {
    let mut prefix = [0_u8; 18];
    let n = file.read(&mut prefix)?;
    #[cfg(feature = "bgzf")]
    if is_bgzf_header(&prefix[..n]) {
        return Ok(InputKind::Bgzf);
    }
    #[cfg(feature = "gzip")]
    if n >= 2 && prefix[..2] == [0x1f, 0x8b] {
        return Ok(InputKind::Gzip);
    }
    Ok(InputKind::Raw)
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
pub fn open_fastq_gzip_libdeflate(path: impl AsRef<Path>) -> Result<FastqReader<Cursor<Vec<u8>>>> {
    open_fastq_gzip_libdeflate_with_config(path, FastqConfig::default())
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
pub fn open_fastq_gzip_libdeflate_with_config(
    path: impl AsRef<Path>,
    config: FastqConfig,
) -> Result<FastqReader<Cursor<Vec<u8>>>> {
    let mut compressed = Vec::new();
    File::open(path)?.read_to_end(&mut compressed)?;
    let decoded = decompress_gzip_libdeflate_buffered(&compressed)?;
    Ok(FastqReader::with_config(Cursor::new(decoded), config))
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
fn decompress_gzip_libdeflate_buffered(compressed: &[u8]) -> Result<Vec<u8>> {
    let mut out = vec![0_u8; initial_gzip_output_capacity(compressed)];
    let mut decompressor = libdeflater::Decompressor::new();
    loop {
        match decompressor.gzip_decompress(compressed, &mut out) {
            Ok(n) => {
                out.truncate(n);
                return Ok(out);
            }
            Err(libdeflater::DecompressionError::InsufficientSpace) => {
                let next = out.len().checked_mul(2).ok_or_else(|| {
                    crate::FastqError::Format("gzip output size exceeds usize range".into())
                })?;
                out.resize(next.max(1), 0);
            }
            Err(err) => {
                return Err(crate::FastqError::Format(format!(
                    "libdeflate gzip inflate failed: {err}"
                )));
            }
        }
    }
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
fn initial_gzip_output_capacity(compressed: &[u8]) -> usize {
    let isize = compressed
        .len()
        .checked_sub(4)
        .and_then(|start| compressed.get(start..))
        .map(|tail| u32::from_le_bytes([tail[0], tail[1], tail[2], tail[3]]) as usize)
        .unwrap_or(0);
    isize.max(compressed.len().saturating_mul(2)).max(1024)
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
pub fn open_fastq_bgzf_parallel_with_options(
    path: impl AsRef<Path>,
    bgzf_config: BgzfParallelConfig,
    fastq_config: FastqConfig,
) -> Result<FastqReader<BgzfParallelReader>> {
    let file = File::open(path)?;
    Ok(FastqReader::with_config(
        BgzfParallelReader::with_config(file, bgzf_config)?,
        fastq_config,
    ))
}

#[cfg(feature = "bgzf")]
pub fn open_fastq_bgzf_adaptive(
    path: impl AsRef<Path>,
    bgzf_config: BgzfParallelConfig,
    fastq_config: FastqConfig,
) -> Result<FastqReader<BgzfAutoReader<File>>> {
    let path = path.as_ref();
    let compressed_len = std::fs::metadata(path)?.len();
    let file = File::open(path)?;
    Ok(FastqReader::with_config(
        BgzfAutoReader::with_config(file, compressed_len, bgzf_config)?,
        fastq_config,
    ))
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

    #[test]
    #[cfg(feature = "libdeflate")]
    fn explicit_libdeflate_gzip_opener_parses_same_file() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "microraptor-libdeflate-gzip-{}.fq.gz",
            std::process::id()
        ));
        let file = File::create(&path).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::fast());
        encoder
            .write_all(b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nJJJJ\n")
            .unwrap();
        encoder.finish().unwrap();

        let mut flate2 = open_fastq(&path).unwrap();
        let mut libdeflate = open_fastq_gzip_libdeflate(&path).unwrap();

        let flate2_stats = crate::benchutil::consume_fastq(&mut flate2).unwrap();
        let libdeflate_stats = crate::benchutil::consume_fastq(&mut libdeflate).unwrap();
        assert_eq!(flate2_stats, libdeflate_stats);

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn paired_openers_parse_matching_files() {
        let dir = std::env::temp_dir();
        let r1_path = dir.join(format!("microraptor-r1-{}.fq", std::process::id()));
        let r2_path = dir.join(format!("microraptor-r2-{}.fq", std::process::id()));
        std::fs::write(&r1_path, b"@frag/1\nACGT\n+\nIIII\n").unwrap();
        std::fs::write(&r2_path, b"@frag/2\nTGCA\n+\nJJJJ\n").unwrap();

        let mut reader = open_paired_fastq(&r1_path, &r2_path).unwrap();
        let batch = reader.next_pair_batch().unwrap().unwrap();
        let pair = batch.pairs().next().unwrap();
        assert_eq!(pair.pair_id(), b"frag");
        assert_eq!(pair.first().seq(), b"ACGT");
        assert_eq!(pair.second().seq(), b"TGCA");

        std::fs::remove_file(r1_path).unwrap();
        std::fs::remove_file(r2_path).unwrap();
    }
}
