use std::fs::File;
#[cfg(all(feature = "gzip", feature = "libdeflate"))]
use std::io::Cursor;
use std::io::Read;
#[cfg(any(feature = "bgzf", feature = "gzip"))]
use std::io::{Seek, SeekFrom};
use std::path::Path;

use crate::error::Result;
use crate::fasta::{FastaConfig, FastaReader};
use crate::fastq::{FastqConfig, FastqReader, PairedFastqReader};
#[cfg(feature = "bgzf")]
use crate::{
    BgzfAutoReader, BgzfInflateBackend, BgzfParallelConfig, BgzfParallelReader, BgzfReader,
    bgzf::is_bgzf_header,
};

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
const DEFAULT_LIBDEFLATE_GZIP_MAX_COMPRESSED_BYTES: usize = if usize::BITS >= 64 {
    1024 * 1024 * 1024
} else {
    usize::MAX / 2
};
#[cfg(all(feature = "gzip", feature = "libdeflate"))]
const DEFAULT_LIBDEFLATE_GZIP_MAX_DECOMPRESSED_BYTES: usize = if usize::BITS >= 64 {
    1024 * 1024 * 1024
} else {
    usize::MAX / 2
};

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
/// Memory limits for explicit buffered libdeflate gzip openers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LibdeflateGzipLimits {
    /// Maximum compressed bytes accepted before buffering.
    pub max_compressed_bytes: usize,
    /// Maximum decompressed bytes accepted before parsing.
    pub max_decompressed_bytes: usize,
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
impl Default for LibdeflateGzipLimits {
    fn default() -> Self {
        Self {
            max_compressed_bytes: DEFAULT_LIBDEFLATE_GZIP_MAX_COMPRESSED_BYTES,
            max_decompressed_bytes: DEFAULT_LIBDEFLATE_GZIP_MAX_DECOMPRESSED_BYTES,
        }
    }
}

#[cfg(any(feature = "bgzf", feature = "gzip"))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum InputKind {
    Raw,
    #[cfg(feature = "gzip")]
    Gzip,
    #[cfg(feature = "bgzf")]
    Bgzf,
}

/// Compression/container detected from input file magic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DetectedInputKind {
    /// Plain uncompressed input.
    Raw,
    /// Ordinary gzip input.
    #[cfg(feature = "gzip")]
    Gzip,
    /// BGZF blocked gzip input.
    #[cfg(feature = "bgzf")]
    Bgzf,
}

#[cfg(any(feature = "bgzf", feature = "gzip"))]
impl From<InputKind> for DetectedInputKind {
    fn from(value: InputKind) -> Self {
        match value {
            InputKind::Raw => Self::Raw,
            #[cfg(feature = "gzip")]
            InputKind::Gzip => Self::Gzip,
            #[cfg(feature = "bgzf")]
            InputKind::Bgzf => Self::Bgzf,
        }
    }
}

/// Detect raw, gzip, or BGZF input by file magic.
pub fn detect_file_input_kind(path: impl AsRef<Path>) -> Result<DetectedInputKind> {
    #[cfg(any(feature = "bgzf", feature = "gzip"))]
    {
        let mut file = File::open(path)?;
        detect_input_kind(&mut file).map(DetectedInputKind::from)
    }
    #[cfg(not(any(feature = "bgzf", feature = "gzip")))]
    {
        let _ = path;
        Ok(DetectedInputKind::Raw)
    }
}

/// Open a FASTQ file with default configuration.
///
/// With default features, this detects raw FASTQ, ordinary gzip, and BGZF by
/// file magic. BGZF is checked before ordinary gzip.
pub fn open_fastq(path: impl AsRef<Path>) -> Result<FastqReader<Box<dyn Read + Send>>> {
    open_fastq_with_config(path, FastqConfig::default())
}

/// Open a FASTQ file with explicit parser configuration.
///
/// This is the primary file-path API for callers that want format
/// auto-detection with custom validation, slab, or pairing settings.
pub fn open_fastq_with_config(
    path: impl AsRef<Path>,
    config: FastqConfig,
) -> Result<FastqReader<Box<dyn Read + Send>>> {
    let reader = open_read_by_magic(path)?;
    Ok(FastqReader::with_config(reader, config))
}

/// Open a FASTA file with default configuration.
///
/// With default features, this detects raw FASTA, ordinary gzip, and BGZF by
/// file magic. BGZF is checked before ordinary gzip.
pub fn open_fasta(path: impl AsRef<Path>) -> Result<FastaReader<Box<dyn Read + Send>>> {
    open_fasta_with_config(path, FastaConfig::default())
}

/// Open a FASTA file with explicit parser configuration.
pub fn open_fasta_with_config(
    path: impl AsRef<Path>,
    config: FastaConfig,
) -> Result<FastaReader<Box<dyn Read + Send>>> {
    let reader = open_read_by_magic(path)?;
    Ok(FastaReader::with_config(reader, config))
}

/// Open a FASTA reference genome with parser settings tuned for long records.
///
/// This uses the same raw/gzip/BGZF transport detection as [`open_fasta`], but
/// lowers records per batch and increases I/O buffering and sequence
/// preallocation hints for chromosome-scale records.
pub fn open_fasta_for_reference(
    path: impl AsRef<Path>,
) -> Result<FastaReader<Box<dyn Read + Send>>> {
    open_fasta_with_config(path, FastaConfig::reference())
}

/// Open ordered R1/R2 FASTQ files with default configuration.
pub fn open_paired_fastq(
    first_path: impl AsRef<Path>,
    second_path: impl AsRef<Path>,
) -> Result<PairedFastqReader<Box<dyn Read + Send>, Box<dyn Read + Send>>> {
    open_paired_fastq_with_config(first_path, second_path, FastqConfig::default())
}

/// Open ordered R1/R2 FASTQ files with the same configuration for both mates.
pub fn open_paired_fastq_with_config(
    first_path: impl AsRef<Path>,
    second_path: impl AsRef<Path>,
    config: FastqConfig,
) -> Result<PairedFastqReader<Box<dyn Read + Send>, Box<dyn Read + Send>>> {
    open_paired_fastq_with_configs(first_path, config.clone(), second_path, config)
}

/// Open ordered R1/R2 FASTQ files with separate mate configurations.
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
            InputKind::Bgzf => {
                let compressed_len = file.metadata()?.len();
                let config = BgzfParallelConfig::new(default_bgzf_workers());
                return Ok(Box::new(BgzfAutoReader::with_config(
                    file,
                    compressed_len,
                    config,
                )?));
            }
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

#[cfg(feature = "bgzf")]
fn default_bgzf_workers() -> usize {
    std::thread::available_parallelism().map_or(1, usize::from)
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
/// Open an ordinary gzip FASTQ file through a buffered libdeflate path.
///
/// Unlike [`open_fastq`], this buffers the fully decompressed input before
/// parsing. It accepts a single gzip member and applies default memory limits;
/// use [`open_fastq_gzip_libdeflate_with_limits`] to set tighter bounds.
pub fn open_fastq_gzip_libdeflate(path: impl AsRef<Path>) -> Result<FastqReader<Cursor<Vec<u8>>>> {
    open_fastq_gzip_libdeflate_with_config(path, FastqConfig::default())
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
/// Open an ordinary gzip FASTQ file through buffered libdeflate with parser
/// configuration.
pub fn open_fastq_gzip_libdeflate_with_config(
    path: impl AsRef<Path>,
    config: FastqConfig,
) -> Result<FastqReader<Cursor<Vec<u8>>>> {
    open_fastq_gzip_libdeflate_with_limits(path, config, LibdeflateGzipLimits::default())
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
/// Open an ordinary single-member gzip FASTQ file through buffered libdeflate
/// with explicit parser configuration and memory limits.
pub fn open_fastq_gzip_libdeflate_with_limits(
    path: impl AsRef<Path>,
    config: FastqConfig,
    limits: LibdeflateGzipLimits,
) -> Result<FastqReader<Cursor<Vec<u8>>>> {
    let compressed = read_limited(path, limits.max_compressed_bytes)?;
    let decoded = decompress_gzip_libdeflate_buffered(&compressed, limits)?;
    Ok(FastqReader::with_config(Cursor::new(decoded), config))
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
/// Open an ordinary gzip FASTA file through a buffered libdeflate path.
///
/// Unlike [`open_fasta`], this buffers the fully decompressed input before
/// parsing. It accepts a single gzip member and applies default memory limits;
/// use [`open_fasta_gzip_libdeflate_with_limits`] to set tighter bounds.
pub fn open_fasta_gzip_libdeflate(path: impl AsRef<Path>) -> Result<FastaReader<Cursor<Vec<u8>>>> {
    open_fasta_gzip_libdeflate_with_config(path, FastaConfig::default())
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
/// Open an ordinary gzip FASTA file through buffered libdeflate with parser
/// configuration.
pub fn open_fasta_gzip_libdeflate_with_config(
    path: impl AsRef<Path>,
    config: FastaConfig,
) -> Result<FastaReader<Cursor<Vec<u8>>>> {
    open_fasta_gzip_libdeflate_with_limits(path, config, LibdeflateGzipLimits::default())
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
/// Open an ordinary single-member gzip FASTA file through buffered libdeflate
/// with explicit parser configuration and memory limits.
pub fn open_fasta_gzip_libdeflate_with_limits(
    path: impl AsRef<Path>,
    config: FastaConfig,
    limits: LibdeflateGzipLimits,
) -> Result<FastaReader<Cursor<Vec<u8>>>> {
    let compressed = read_limited(path, limits.max_compressed_bytes)?;
    let decoded = decompress_gzip_libdeflate_buffered(&compressed, limits)?;
    Ok(FastaReader::with_config(Cursor::new(decoded), config))
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
fn read_limited(path: impl AsRef<Path>, max_bytes: usize) -> Result<Vec<u8>> {
    let path = path.as_ref();
    let len = std::fs::metadata(path)?.len();
    if len > max_bytes as u64 {
        return Err(crate::FastqError::Format(format!(
            "gzip input exceeds libdeflate compressed limit ({len} > {max_bytes} bytes)"
        )));
    }
    let mut compressed = Vec::with_capacity(usize::try_from(len).unwrap_or(0));
    File::open(path)?.read_to_end(&mut compressed)?;
    if compressed.len() > max_bytes {
        return Err(crate::FastqError::Format(format!(
            "gzip input exceeds libdeflate compressed limit ({} > {max_bytes} bytes)",
            compressed.len()
        )));
    }
    Ok(compressed)
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
fn decompress_gzip_libdeflate_buffered(
    compressed: &[u8],
    limits: LibdeflateGzipLimits,
) -> Result<Vec<u8>> {
    if limits.max_decompressed_bytes == 0 {
        return Err(crate::FastqError::Format(
            "libdeflate decompressed limit must be greater than zero".into(),
        ));
    }
    ensure_single_gzip_member(compressed)?;
    let initial = initial_gzip_output_capacity(compressed).min(limits.max_decompressed_bytes);
    let mut out = vec![0_u8; initial.max(1)];
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
                if next > limits.max_decompressed_bytes {
                    return Err(crate::FastqError::Format(format!(
                        "gzip output exceeds libdeflate decompressed limit (>{} bytes)",
                        limits.max_decompressed_bytes
                    )));
                }
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
fn ensure_single_gzip_member(compressed: &[u8]) -> Result<()> {
    let end = gzip_first_member_end(compressed)?;
    if end != compressed.len() {
        return Err(crate::FastqError::Format(
            "libdeflate gzip opener accepts exactly one gzip member; use open_fastq/open_fasta for concatenated gzip streams".into(),
        ));
    }
    Ok(())
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
fn gzip_first_member_end(compressed: &[u8]) -> Result<usize> {
    use flate2::{Decompress, FlushDecompress, Status};

    if compressed.len() < 18 || compressed[..2] != [0x1f, 0x8b] || compressed[2] != 8 {
        return Err(crate::FastqError::Format("invalid gzip header".into()));
    }
    let flags = compressed[3];
    if flags & 0xe0 != 0 {
        return Err(crate::FastqError::Format(
            "gzip header uses reserved flags".into(),
        ));
    }

    let mut pos = 10;
    if flags & 0x04 != 0 {
        let xlen = read_gzip_u16(compressed, pos)? as usize;
        pos = pos
            .checked_add(2)
            .and_then(|p| p.checked_add(xlen))
            .ok_or_else(|| crate::FastqError::Format("gzip extra field is too large".into()))?;
        if pos > compressed.len() {
            return Err(crate::FastqError::Format(
                "truncated gzip extra field".into(),
            ));
        }
    }
    if flags & 0x08 != 0 {
        pos = scan_gzip_cstring(compressed, pos, "name")?;
    }
    if flags & 0x10 != 0 {
        pos = scan_gzip_cstring(compressed, pos, "comment")?;
    }
    if flags & 0x02 != 0 {
        pos = pos
            .checked_add(2)
            .ok_or_else(|| crate::FastqError::Format("gzip header CRC is too large".into()))?;
        if pos > compressed.len() {
            return Err(crate::FastqError::Format(
                "truncated gzip header CRC".into(),
            ));
        }
    }
    if pos + 8 > compressed.len() {
        return Err(crate::FastqError::Format("truncated gzip member".into()));
    }

    let deflate = &compressed[pos..];
    let mut decoder = Decompress::new(false);
    let mut scratch = [0_u8; 8192];
    loop {
        let consumed_before = decoder.total_in();
        let produced_before = decoder.total_out();
        let consumed = usize::try_from(consumed_before)
            .map_err(|_| crate::FastqError::Format("gzip member exceeds usize range".into()))?;
        let status = decoder
            .decompress(
                deflate.get(consumed..).ok_or_else(|| {
                    crate::FastqError::Format("truncated gzip deflate stream".into())
                })?,
                &mut scratch,
                FlushDecompress::None,
            )
            .map_err(|err| {
                crate::FastqError::Format(format!("gzip deflate parse failed: {err}"))
            })?;
        if status == Status::StreamEnd {
            let deflate_len = usize::try_from(decoder.total_in()).map_err(|_| {
                crate::FastqError::Format("gzip deflate stream exceeds usize range".into())
            })?;
            let end = pos
                .checked_add(deflate_len)
                .and_then(|p| p.checked_add(8))
                .ok_or_else(|| crate::FastqError::Format("gzip member is too large".into()))?;
            if end > compressed.len() {
                return Err(crate::FastqError::Format("truncated gzip trailer".into()));
            }
            return Ok(end);
        }
        if decoder.total_in() == consumed_before && decoder.total_out() == produced_before {
            return Err(crate::FastqError::Format(
                "gzip deflate parser made no progress".into(),
            ));
        }
    }
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
fn read_gzip_u16(bytes: &[u8], pos: usize) -> Result<u16> {
    let Some(pair) = bytes.get(pos..pos + 2) else {
        return Err(crate::FastqError::Format("truncated gzip header".into()));
    };
    Ok(u16::from_le_bytes([pair[0], pair[1]]))
}

#[cfg(all(feature = "gzip", feature = "libdeflate"))]
fn scan_gzip_cstring(bytes: &[u8], pos: usize, field: &str) -> Result<usize> {
    let Some(relative_end) = bytes[pos..].iter().position(|&b| b == 0) else {
        return Err(crate::FastqError::Format(format!(
            "unterminated gzip {field} field"
        )));
    };
    Ok(pos + relative_end + 1)
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
/// Open a BGZF FASTQ file with an explicit parallel worker count.
pub fn open_fastq_bgzf_parallel(
    path: impl AsRef<Path>,
    workers: usize,
) -> Result<FastqReader<BgzfParallelReader>> {
    open_fastq_bgzf_parallel_with_config(path, workers, FastqConfig::default())
}

#[cfg(feature = "bgzf")]
/// Open a BGZF FASTQ file with an explicit parallel worker count and parser
/// configuration.
pub fn open_fastq_bgzf_parallel_with_config(
    path: impl AsRef<Path>,
    workers: usize,
    config: FastqConfig,
) -> Result<FastqReader<BgzfParallelReader>> {
    open_fastq_bgzf_parallel_with_backend(path, workers, BgzfInflateBackend::default(), config)
}

#[cfg(feature = "bgzf")]
/// Open a BGZF FASTQ file with explicit BGZF and FASTQ configuration.
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
/// Open a BGZF FASTQ file with adaptive serial/parallel reading.
///
/// The supplied [`BgzfParallelConfig`] controls the parallelization threshold,
/// backend, workers, queue depths, and optional metrics.
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
/// Open a BGZF FASTQ file through the flate2 serial backend.
pub fn open_fastq_bgzf_flate2(path: impl AsRef<Path>) -> Result<FastqReader<BgzfReader<File>>> {
    open_fastq_bgzf_with_backend(path, BgzfInflateBackend::Flate2, FastqConfig::default())
}

#[cfg(feature = "bgzf")]
/// Open a BGZF FASTQ file with a selected serial inflate backend.
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
/// Open a BGZF FASTQ file with a selected parallel inflate backend.
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

    #[cfg(all(feature = "gzip", feature = "libdeflate"))]
    fn gzip_member(payload: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        {
            let mut encoder = flate2::write::GzEncoder::new(&mut out, flate2::Compression::fast());
            encoder.write_all(payload).unwrap();
            encoder.finish().unwrap();
        }
        out
    }

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
    #[cfg(feature = "gzip")]
    fn opens_gzip_fasta_by_magic() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("microraptor-{}.fa.gz", std::process::id()));
        let file = File::create(&path).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::fast());
        encoder.write_all(b">seq1\nAC\nGT\n").unwrap();
        encoder.finish().unwrap();

        let mut reader = open_fasta(&path).unwrap();
        let batch = reader.next_batch().unwrap().unwrap();
        let rec = batch.records().next().unwrap();
        assert_eq!(rec.name_without_gt(), b"seq1");
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
    #[cfg(feature = "bgzf")]
    fn opens_bgzf_fasta_by_magic() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("microraptor-fasta-{}.bgz", std::process::id()));
        let mut writer = crate::BgzfWriter::new(File::create(&path).unwrap());
        writer.write_all(b">seq1\nAC\nGT\n").unwrap();
        writer.finish().unwrap();

        let mut reader = open_fasta(&path).unwrap();
        let batch = reader.next_batch().unwrap().unwrap();
        let rec = batch.records().next().unwrap();
        assert_eq!(rec.name_without_gt(), b"seq1");
        assert_eq!(rec.seq(), b"ACGT");

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    #[cfg(feature = "bgzf")]
    fn public_bgzf_adaptive_opener_selects_serial_and_parallel() {
        let dir = std::env::temp_dir();
        let small_path = dir.join(format!("microraptor-bgzf-small-{}.bgz", std::process::id()));
        let large_path = dir.join(format!("microraptor-bgzf-large-{}.bgz", std::process::id()));

        let mut small_writer = crate::BgzfWriter::new(File::create(&small_path).unwrap());
        small_writer.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
        small_writer.finish().unwrap();

        let mut large_writer = crate::BgzfWriter::new(File::create(&large_path).unwrap());
        for i in 0..2048 {
            writeln!(large_writer, "@r{i}\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII").unwrap();
        }
        large_writer.finish().unwrap();

        let serial = open_fastq_bgzf_adaptive(
            &small_path,
            BgzfParallelConfig::new(2).with_parallel_min_compressed_bytes(u64::MAX),
            FastqConfig::default(),
        )
        .unwrap()
        .into_inner();
        assert!(matches!(serial, crate::BgzfAutoReader::Serial(_)));

        let parallel = open_fastq_bgzf_adaptive(
            &large_path,
            BgzfParallelConfig::new(2).with_parallel_min_compressed_bytes(0),
            FastqConfig::default(),
        )
        .unwrap()
        .into_inner();
        assert!(matches!(parallel, crate::BgzfAutoReader::Parallel(_)));

        std::fs::remove_file(small_path).unwrap();
        std::fs::remove_file(large_path).unwrap();
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
    #[cfg(feature = "libdeflate")]
    fn explicit_libdeflate_gzip_rejects_concatenated_fastq_members() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "microraptor-libdeflate-gzip-multi-{}.fq.gz",
            std::process::id()
        ));
        let mut encoded = gzip_member(b"@r1\nACGT\n+\nIIII\n");
        encoded.extend_from_slice(&gzip_member(b"@r2\nTGCA\n+\nJJJJ\n"));
        std::fs::write(&path, encoded).unwrap();

        let mut default_reader = open_fastq(&path).unwrap();
        let default_stats = crate::benchutil::consume_fastq(&mut default_reader).unwrap();
        assert_eq!(default_stats.records, 2);

        let err = match open_fastq_gzip_libdeflate(&path) {
            Ok(_) => panic!("concatenated gzip member was accepted"),
            Err(err) => err,
        };
        assert!(err.to_string().contains("one gzip member"));

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    #[cfg(feature = "libdeflate")]
    fn explicit_libdeflate_gzip_limits_compressed_input() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "microraptor-libdeflate-gzip-limit-{}.fq.gz",
            std::process::id()
        ));
        std::fs::write(&path, gzip_member(b"@r1\nACGT\n+\nIIII\n")).unwrap();

        let err = match open_fastq_gzip_libdeflate_with_limits(
            &path,
            FastqConfig::default(),
            LibdeflateGzipLimits {
                max_compressed_bytes: 1,
                max_decompressed_bytes: 1024,
            },
        ) {
            Ok(_) => panic!("compressed limit was not enforced"),
            Err(err) => err,
        };
        assert!(err.to_string().contains("compressed limit"));

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    #[cfg(feature = "libdeflate")]
    fn explicit_libdeflate_gzip_limits_decompressed_output() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "microraptor-libdeflate-gzip-output-limit-{}.fq.gz",
            std::process::id()
        ));
        std::fs::write(&path, gzip_member(b"@r1\nACGT\n+\nIIII\n")).unwrap();

        let err = match open_fastq_gzip_libdeflate_with_limits(
            &path,
            FastqConfig::default(),
            LibdeflateGzipLimits {
                max_compressed_bytes: 1024,
                max_decompressed_bytes: 4,
            },
        ) {
            Ok(_) => panic!("decompressed limit was not enforced"),
            Err(err) => err,
        };
        assert!(err.to_string().contains("decompressed limit"));

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    #[cfg(feature = "libdeflate")]
    fn explicit_libdeflate_gzip_fasta_opener_parses_same_file() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "microraptor-libdeflate-gzip-fasta-{}.fa.gz",
            std::process::id()
        ));
        let file = File::create(&path).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::fast());
        encoder.write_all(b">seq1\nAC\nGT\n").unwrap();
        encoder.finish().unwrap();

        let mut flate2 = open_fasta(&path).unwrap();
        let mut libdeflate = open_fasta_gzip_libdeflate(&path).unwrap();

        let flate2_stats = crate::benchutil::consume_fasta(&mut flate2).unwrap();
        let libdeflate_stats = crate::benchutil::consume_fasta(&mut libdeflate).unwrap();
        assert_eq!(flate2_stats, libdeflate_stats);

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    #[cfg(feature = "libdeflate")]
    fn explicit_libdeflate_gzip_rejects_concatenated_fasta_members() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "microraptor-libdeflate-gzip-fasta-multi-{}.fa.gz",
            std::process::id()
        ));
        let mut encoded = gzip_member(b">seq1\nAC\n");
        encoded.extend_from_slice(&gzip_member(b">seq2\nGT\n"));
        std::fs::write(&path, encoded).unwrap();

        let mut default_reader = open_fasta(&path).unwrap();
        let default_stats = crate::benchutil::consume_fasta(&mut default_reader).unwrap();
        assert_eq!(default_stats.records, 2);

        let err = match open_fasta_gzip_libdeflate(&path) {
            Ok(_) => panic!("concatenated gzip member was accepted"),
            Err(err) => err,
        };
        assert!(err.to_string().contains("one gzip member"));

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

    #[test]
    fn open_fasta_parses_raw_file() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("microraptor-{}.fa", std::process::id()));
        std::fs::write(&path, b">seq1\nAC\nGT\n").unwrap();

        let mut reader = open_fasta(&path).unwrap();
        let batch = reader.next_batch().unwrap().unwrap();
        let rec = batch.records().next().unwrap();
        assert_eq!(rec.id_token(), b"seq1");
        assert_eq!(rec.seq(), b"ACGT");

        std::fs::remove_file(path).unwrap();
    }
}
