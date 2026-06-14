#[cfg(feature = "gzip")]
use std::io::Write;
use std::path::{Path, PathBuf};
#[cfg(feature = "bgzf")]
use std::sync::Arc;
use std::time::{Duration, Instant};

use microraptor::benchutil::{
    StreamStats, consume_fasta, consume_fastq, consume_trusted_fastq_read_direct_with_pack,
    consume_trusted_fastq_read_with_pack, synthetic_fasta, synthetic_fastq,
};
use microraptor::pack::{TrustedPackedRecord, pack_bases_and_summarize_qualities_into};
use microraptor::{
    DetectedInputKind, FastaReader, FastqConfig, FastqReader, PairValidation, Result,
};

#[path = "microraptor_bench/report.rs"]
mod report;

enum BenchRead {
    Raw(std::fs::File),
    #[cfg(feature = "gzip")]
    Gzip(flate2::read::MultiGzDecoder<std::fs::File>),
    #[cfg(feature = "bgzf")]
    Bgzf(microraptor::BgzfAutoReader<std::fs::File>),
}

impl std::io::Read for BenchRead {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        match self {
            Self::Raw(reader) => reader.read(out),
            #[cfg(feature = "gzip")]
            Self::Gzip(reader) => reader.read(out),
            #[cfg(feature = "bgzf")]
            Self::Bgzf(reader) => reader.read(out),
        }
    }
}

#[derive(Debug, Clone)]
struct Config {
    records: usize,
    read_len: usize,
    iters: usize,
    slab_size: usize,
    workers: usize,
    bgzf_parallel_min_bytes: u64,
    json: bool,
    input: Option<PathBuf>,
    paired_inputs: Option<(PathBuf, PathBuf)>,
    mode: Mode,
    format: InputFormat,
    bgzf_pack_check: Option<BgzfPackCheck>,
    profile_bgzf_parallel: bool,
    mmap: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            records: 200_000,
            read_len: 150,
            iters: 5,
            slab_size: 256 * 1024,
            workers: std::thread::available_parallelism().map_or(1, usize::from),
            bgzf_parallel_min_bytes: 512 * 1024,
            json: false,
            input: None,
            paired_inputs: None,
            mode: Mode::All,
            format: InputFormat::Fastq,
            bgzf_pack_check: None,
            profile_bgzf_parallel: false,
            mmap: false,
        }
    }
}

#[derive(Debug, Clone)]
struct BgzfPackCheck {
    label: String,
    min_input_bytes: usize,
    tolerance_pct: u128,
    check_timing: bool,
}

impl Default for BgzfPackCheck {
    fn default() -> Self {
        Self {
            label: "bgzf".into(),
            min_input_bytes: 0,
            tolerance_pct: 15,
            check_timing: true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    All,
    Parse,
    Pack,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum InputFormat {
    Fastq,
    Fasta,
}

impl InputFormat {
    fn as_str(self) -> &'static str {
        match self {
            Self::Fastq => "fastq",
            Self::Fasta => "fasta",
        }
    }
}

impl Mode {
    fn includes_parse(self) -> bool {
        matches!(self, Self::All | Self::Parse)
    }

    fn includes_pack(self) -> bool {
        matches!(self, Self::All | Self::Pack)
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::Parse => "parse",
            Self::Pack => "pack",
        }
    }
}

#[derive(Debug, Clone)]
struct Measurement {
    name: String,
    bytes: usize,
    records: u64,
    bases: u64,
    best: Duration,
    samples: Vec<Duration>,
    checksum: u64,
    extras: Vec<(&'static str, u64)>,
}

fn main() {
    match run() {
        Ok(()) => {}
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    }
}

fn run() -> Result<()> {
    let config = parse_args();
    if let Some((first, second)) = config.paired_inputs.as_ref() {
        return run_paired_input(first, second, &config);
    }
    if let Some(path) = config.input.as_deref() {
        return run_real_input(path, &config);
    }

    if config.format == InputFormat::Fasta {
        return run_synthetic_fasta(&config);
    }

    let raw = synthetic_fastq(config.records, config.read_len);
    #[cfg_attr(not(any(feature = "gzip", feature = "bgzf")), allow(unused_mut))]
    let mut measurements = Vec::new();
    if config.mode.includes_parse() {
        measurements.push(measure_fastq("raw", &raw, &config)?);
    }
    if config.mode.includes_pack() {
        measurements.push(measure_trusted_pack("pack-seq-qual", &raw, &config)?);
        measurements.push(measure_direct_pack("direct-pack-seq-qual", &raw, &config)?);
        measurements.push(measure_pack("reader-pack-seq-qual", &raw, &config)?);
    }

    #[cfg(feature = "gzip")]
    {
        if config.mode.includes_parse() {
            let gzip = gzip_bytes(&raw)?;
            measurements.push(measure_gzip("gzip", &gzip, &config)?);
        }
    }

    #[cfg(feature = "bgzf")]
    {
        if config.mode.includes_parse() {
            let bgzf = microraptor::compress_bgzf_parallel(&raw, config.workers)?;
            measurements.push(measure_bgzf_serial("bgzf-serial", &bgzf, &config)?);
            measurements.push(measure_bgzf_parallel("bgzf-parallel", &bgzf, &config)?);
            #[cfg(feature = "libdeflate")]
            {
                measurements.push(measure_bgzf_libdeflate_serial(
                    "bgzf-libdeflate-serial",
                    &bgzf,
                    &config,
                )?);
                measurements.push(measure_bgzf_libdeflate_parallel(
                    "bgzf-libdeflate-parallel",
                    &bgzf,
                    &config,
                )?);
            }
        }
        if config.mode.includes_pack() {
            let bgzf = microraptor::compress_bgzf_parallel(&raw, config.workers)?;
            measurements.push(measure_bgzf_trusted_pack(
                "bgzf-pack-seq-qual",
                &bgzf,
                &config,
            )?);
            measurements.push(measure_bgzf_adaptive_trusted_pack(
                "bgzf-adaptive-pack-seq-qual",
                &bgzf,
                &config,
            )?);
            #[cfg(feature = "libdeflate")]
            {
                measurements.push(measure_bgzf_libdeflate_serial_trusted_pack(
                    "bgzf-libdeflate-serial-pack-seq-qual",
                    &bgzf,
                    &config,
                )?);
                measurements.push(measure_bgzf_libdeflate_adaptive_trusted_pack(
                    "bgzf-libdeflate-adaptive-pack-seq-qual",
                    &bgzf,
                    &config,
                )?);
            }
        }
    }

    emit_report(&config, "synthetic", raw.len(), &measurements);
    Ok(())
}

fn run_synthetic_fasta(config: &Config) -> Result<()> {
    let raw = synthetic_fasta(config.records, config.read_len);
    #[cfg_attr(not(any(feature = "gzip", feature = "bgzf")), allow(unused_mut))]
    let mut measurements = vec![measure_fasta("fasta-raw", &raw, config)?];

    #[cfg(feature = "gzip")]
    {
        let gzip = gzip_bytes(&raw)?;
        measurements.push(measure_gzip_fasta("fasta-gzip", &gzip, config)?);
    }

    #[cfg(feature = "bgzf")]
    {
        let bgzf = microraptor::compress_bgzf_parallel(&raw, config.workers)?;
        measurements.push(measure_bgzf_fasta_serial(
            "fasta-bgzf-serial",
            &bgzf,
            config,
        )?);
        measurements.push(measure_bgzf_fasta_parallel(
            "fasta-bgzf-parallel",
            &bgzf,
            config,
        )?);
    }

    emit_report(config, "synthetic-fasta", raw.len(), &measurements);
    Ok(())
}

fn run_paired_input(first: &Path, second: &Path, config: &Config) -> Result<()> {
    let input_bytes = checked_file_len(first)?
        .checked_add(checked_file_len(second)?)
        .ok_or_else(|| microraptor::FastqError::Format("paired input size overflow".into()))?;
    let mut measurements = Vec::new();
    if config.mode.includes_parse() {
        measurements.push(measure_paired_path_fastq(
            "file-paired-auto",
            first,
            second,
            input_bytes,
            config,
        )?);
    }
    if config.mode.includes_pack() {
        measurements.push(measure_paired_path_pack(
            "file-paired-pack-seq-qual",
            first,
            second,
            input_bytes,
            config,
        )?);
        measurements.push(measure_paired_path_trusted_pack(
            "file-paired-trusted-pack-seq-qual",
            first,
            second,
            input_bytes,
            config,
        )?);
    }

    let source = format!("{}+{}", first.display(), second.display());
    emit_report(config, &source, input_bytes, &measurements);
    Ok(())
}

fn run_real_input(path: &Path, config: &Config) -> Result<()> {
    let input_bytes = checked_file_len(path)?;
    let mut measurements = Vec::new();
    if config.format == InputFormat::Fasta {
        measurements.push(measure_path_fasta(
            "file-fasta-auto",
            path,
            input_bytes,
            config,
        )?);
        #[cfg(feature = "mmap")]
        if config.mmap && microraptor::detect_file_input_kind(path)? == DetectedInputKind::Raw {
            measurements.push(measure_path_fasta_mmap(
                "file-fasta-mmap",
                path,
                input_bytes,
                config,
            )?);
        }
        let source = path.to_string_lossy();
        emit_report(config, &source, input_bytes, &measurements);
        return Ok(());
    }

    if config.mode.includes_parse() {
        measurements.push(measure_path_fastq("file-auto", path, input_bytes, config)?);
        #[cfg(feature = "mmap")]
        if config.mmap && microraptor::detect_file_input_kind(path)? == DetectedInputKind::Raw {
            measurements.push(measure_path_fastq_mmap(
                "file-mmap",
                path,
                input_bytes,
                config,
            )?);
        }
        #[cfg(all(feature = "bgzf", feature = "libdeflate"))]
        if path_has_bgzf_header(path)? {
            measurements.push(measure_path_bgzf_libdeflate_serial(
                "file-bgzf-libdeflate-serial",
                path,
                input_bytes,
                config,
            )?);
            measurements.push(measure_path_bgzf_libdeflate_parallel(
                "file-bgzf-libdeflate-parallel",
                path,
                input_bytes,
                config,
            )?);
        }
    }
    if config.mode.includes_pack() {
        measurements.push(measure_path_trusted_pack(
            "file-pack-seq-qual",
            path,
            input_bytes,
            config,
        )?);
        measurements.push(measure_path_direct_pack(
            "file-direct-pack-seq-qual",
            path,
            input_bytes,
            config,
        )?);
        measurements.push(measure_path_pack(
            "file-reader-pack-seq-qual",
            path,
            input_bytes,
            config,
        )?);
        #[cfg(feature = "bgzf")]
        if path_has_bgzf_header(path)? {
            measurements.push(measure_path_bgzf_adaptive_trusted_pack(
                "file-bgzf-adaptive-pack-seq-qual",
                path,
                input_bytes,
                config,
            )?);
            #[cfg(feature = "libdeflate")]
            {
                measurements.push(measure_path_bgzf_libdeflate_serial_trusted_pack(
                    "file-bgzf-libdeflate-serial-pack-seq-qual",
                    path,
                    input_bytes,
                    config,
                )?);
                measurements.push(measure_path_bgzf_libdeflate_adaptive_trusted_pack(
                    "file-bgzf-libdeflate-adaptive-pack-seq-qual",
                    path,
                    input_bytes,
                    config,
                )?);
            }
            if config.profile_bgzf_parallel {
                measurements.push(measure_path_bgzf_direct_parallel_trusted_pack(
                    "file-bgzf-direct-parallel-pack-seq-qual",
                    path,
                    input_bytes,
                    config,
                )?);
            }
        }
    }

    if let Some(check) = config.bgzf_pack_check.as_ref() {
        return check_bgzf_pack_regression(check, &measurements);
    }

    let source = path.to_string_lossy();

    emit_report(config, &source, input_bytes, &measurements);
    Ok(())
}

fn check_bgzf_pack_regression(check: &BgzfPackCheck, rows: &[Measurement]) -> Result<()> {
    let reader = required_measurement(rows, "file-reader-pack-seq-qual")?;
    let default = required_measurement(rows, "file-pack-seq-qual")?;
    let adaptive = required_measurement(rows, "file-bgzf-adaptive-pack-seq-qual")?;
    let libdeflate_serial =
        required_measurement(rows, "file-bgzf-libdeflate-serial-pack-seq-qual")?;
    let libdeflate_adaptive =
        required_measurement(rows, "file-bgzf-libdeflate-adaptive-pack-seq-qual")?;

    if default.bytes < check.min_input_bytes {
        return Err(microraptor::FastqError::Format(format!(
            "BGZF {} fixture is below parallel threshold: {} bytes < {} bytes",
            check.label, default.bytes, check.min_input_bytes
        )));
    }

    for row in [default, adaptive, libdeflate_serial, libdeflate_adaptive] {
        if row.checksum != reader.checksum {
            return Err(microraptor::FastqError::Format(format!(
                "BGZF {} checksum mismatch: {}={} reader={}",
                check.label, row.name, row.checksum, reader.checksum
            )));
        }
    }

    if check.check_timing {
        let reader_limit = tolerance_limit(reader.best.as_nanos(), check.tolerance_pct);
        if adaptive.best.as_nanos() > reader_limit {
            return Err(microraptor::FastqError::Format(format!(
                "BGZF {} adaptive pack regression: {} ns > {} ns",
                check.label,
                adaptive.best.as_nanos(),
                reader_limit
            )));
        }

        let adaptive_limit = tolerance_limit(adaptive.best.as_nanos(), check.tolerance_pct);
        if default.best.as_nanos() > adaptive_limit {
            return Err(microraptor::FastqError::Format(format!(
                "BGZF {} default pack is slower than adaptive: {} ns > {} ns",
                check.label,
                default.best.as_nanos(),
                adaptive_limit
            )));
        }

        if libdeflate_adaptive.best.as_nanos() > reader_limit {
            return Err(microraptor::FastqError::Format(format!(
                "BGZF {} libdeflate adaptive pack regression: {} ns > {} ns",
                check.label,
                libdeflate_adaptive.best.as_nanos(),
                reader_limit
            )));
        }
    }

    println!("{}_input_bytes\t{}", check.label, default.bytes);
    println!(
        "{}_file-pack-seq-qual_ns\t{}",
        check.label,
        default.best.as_nanos()
    );
    println!(
        "{}_file-reader-pack-seq-qual_ns\t{}",
        check.label,
        reader.best.as_nanos()
    );
    println!(
        "{}_file-bgzf-adaptive-pack-seq-qual_ns\t{}",
        check.label,
        adaptive.best.as_nanos()
    );
    println!(
        "{}_file-bgzf-libdeflate-serial-pack-seq-qual_ns\t{}",
        check.label,
        libdeflate_serial.best.as_nanos()
    );
    println!(
        "{}_file-bgzf-libdeflate-adaptive-pack-seq-qual_ns\t{}",
        check.label,
        libdeflate_adaptive.best.as_nanos()
    );
    println!("{}_timing_checks\t{}", check.label, check.check_timing);
    Ok(())
}

fn required_measurement<'a>(rows: &'a [Measurement], name: &str) -> Result<&'a Measurement> {
    rows.iter().find(|row| row.name == name).ok_or_else(|| {
        microraptor::FastqError::Format(format!("missing BGZF pack benchmark row: {name}"))
    })
}

fn tolerance_limit(ns: u128, tolerance_pct: u128) -> u128 {
    ns + (ns * tolerance_pct / 100)
}

fn checked_file_len(path: &Path) -> Result<usize> {
    usize::try_from(std::fs::metadata(path)?.len())
        .map_err(|_| microraptor::FastqError::Format("input file is too large".into()))
}

fn open_bench_read(path: &Path, _config: &Config) -> Result<BenchRead> {
    match microraptor::detect_file_input_kind(path)? {
        #[cfg(feature = "bgzf")]
        DetectedInputKind::Bgzf => {
            let file = std::fs::File::open(path)?;
            let compressed_len = file.metadata()?.len();
            Ok(BenchRead::Bgzf(microraptor::BgzfAutoReader::with_config(
                file,
                compressed_len,
                bgzf_config(_config),
            )?))
        }
        #[cfg(feature = "gzip")]
        DetectedInputKind::Gzip => {
            let file = std::fs::File::open(path)?;
            Ok(BenchRead::Gzip(flate2::read::MultiGzDecoder::new(file)))
        }
        DetectedInputKind::Raw => Ok(BenchRead::Raw(std::fs::File::open(path)?)),
    }
}

fn measure_fastq(name: &str, input: &[u8], config: &Config) -> Result<Measurement> {
    measure(name, input.len(), config.iters, || {
        let source = std::io::Cursor::new(input);
        let mut reader = FastqReader::with_config(
            source,
            FastqConfig {
                slab_size: config.slab_size,
                validate: true,
                ..FastqConfig::default()
            },
        );
        consume_fastq(&mut reader)
    })
}

fn measure_fasta(name: &str, input: &[u8], config: &Config) -> Result<Measurement> {
    measure(name, input.len(), config.iters, || {
        let source = std::io::Cursor::new(input);
        let mut reader = FastaReader::with_config(
            source,
            microraptor::FastaConfig {
                batch_records: fasta_batch_records(config),
                ..microraptor::FastaConfig::default()
            },
        );
        consume_fasta(&mut reader)
    })
}

#[cfg(feature = "gzip")]
fn measure_gzip(name: &str, input: &[u8], config: &Config) -> Result<Measurement> {
    measure(name, input.len(), config.iters, || {
        let source = flate2::read::MultiGzDecoder::new(input);
        let mut reader = FastqReader::with_config(
            source,
            FastqConfig {
                slab_size: config.slab_size,
                validate: true,
                ..FastqConfig::default()
            },
        );
        consume_fastq(&mut reader)
    })
}

#[cfg(feature = "gzip")]
fn measure_gzip_fasta(name: &str, input: &[u8], config: &Config) -> Result<Measurement> {
    measure(name, input.len(), config.iters, || {
        let source = flate2::read::MultiGzDecoder::new(input);
        let mut reader = FastaReader::with_config(
            source,
            microraptor::FastaConfig {
                batch_records: fasta_batch_records(config),
                ..microraptor::FastaConfig::default()
            },
        );
        consume_fasta(&mut reader)
    })
}

#[cfg(feature = "bgzf")]
fn measure_bgzf_serial(name: &str, input: &[u8], config: &Config) -> Result<Measurement> {
    measure(name, input.len(), config.iters, || {
        let source = microraptor::BgzfReader::new(input);
        let mut reader = FastqReader::with_config(
            source,
            FastqConfig {
                slab_size: config.slab_size,
                validate: true,
                ..FastqConfig::default()
            },
        );
        consume_fastq(&mut reader)
    })
}

#[cfg(feature = "bgzf")]
fn measure_bgzf_fasta_serial(name: &str, input: &[u8], config: &Config) -> Result<Measurement> {
    measure(name, input.len(), config.iters, || {
        let source = microraptor::BgzfReader::new(input);
        let mut reader = FastaReader::with_config(
            source,
            microraptor::FastaConfig {
                batch_records: fasta_batch_records(config),
                ..microraptor::FastaConfig::default()
            },
        );
        consume_fasta(&mut reader)
    })
}

#[cfg(feature = "bgzf")]
fn measure_bgzf_parallel(name: &str, input: &[u8], config: &Config) -> Result<Measurement> {
    let owned: Arc<[u8]> = Arc::from(input);
    measure(name, input.len(), config.iters, || {
        let source = microraptor::BgzfParallelReader::new(
            std::io::Cursor::new(Arc::clone(&owned)),
            config.workers,
        )?;
        let mut reader = FastqReader::with_config(
            source,
            FastqConfig {
                slab_size: config.slab_size,
                validate: true,
                ..FastqConfig::default()
            },
        );
        consume_fastq(&mut reader)
    })
}

#[cfg(feature = "bgzf")]
fn measure_bgzf_fasta_parallel(name: &str, input: &[u8], config: &Config) -> Result<Measurement> {
    let owned: Arc<[u8]> = Arc::from(input);
    measure(name, input.len(), config.iters, || {
        let source = microraptor::BgzfParallelReader::new(
            std::io::Cursor::new(Arc::clone(&owned)),
            config.workers,
        )?;
        let mut reader = FastaReader::with_config(
            source,
            microraptor::FastaConfig {
                batch_records: fasta_batch_records(config),
                ..microraptor::FastaConfig::default()
            },
        );
        consume_fasta(&mut reader)
    })
}

#[cfg(feature = "bgzf")]
fn measure_bgzf_trusted_pack(name: &str, input: &[u8], config: &Config) -> Result<Measurement> {
    measure(name, input.len(), config.iters, || {
        let source = microraptor::BgzfReader::new(input);
        consume_trusted_fastq_read_with_pack(source, fastq_config(config))
    })
}

#[cfg(feature = "bgzf")]
fn measure_bgzf_adaptive_trusted_pack(
    name: &str,
    input: &[u8],
    config: &Config,
) -> Result<Measurement> {
    let owned: Arc<[u8]> = Arc::from(input);
    measure_with_bgzf_metrics(name, input.len(), config.iters, |metrics| {
        let source = microraptor::BgzfAutoReader::with_config(
            std::io::Cursor::new(Arc::clone(&owned)),
            input.len() as u64,
            bgzf_config(config).with_metrics(metrics),
        )?;
        consume_trusted_fastq_read_with_pack(source, fastq_config(config))
    })
}

#[cfg(all(feature = "bgzf", feature = "libdeflate"))]
fn measure_bgzf_libdeflate_serial(
    name: &str,
    input: &[u8],
    config: &Config,
) -> Result<Measurement> {
    measure(name, input.len(), config.iters, || {
        let source = microraptor::BgzfReader::with_inflate_backend(
            input,
            microraptor::BgzfInflateBackend::Libdeflate,
        );
        let mut reader = FastqReader::with_config(
            source,
            FastqConfig {
                slab_size: config.slab_size,
                validate: true,
                ..FastqConfig::default()
            },
        );
        consume_fastq(&mut reader)
    })
}

#[cfg(all(feature = "bgzf", feature = "libdeflate"))]
fn measure_bgzf_libdeflate_parallel(
    name: &str,
    input: &[u8],
    config: &Config,
) -> Result<Measurement> {
    let owned: Arc<[u8]> = Arc::from(input);
    measure(name, input.len(), config.iters, || {
        let source = microraptor::BgzfParallelReader::with_inflate_backend(
            std::io::Cursor::new(Arc::clone(&owned)),
            config.workers,
            microraptor::BgzfInflateBackend::Libdeflate,
        )?;
        let mut reader = FastqReader::with_config(
            source,
            FastqConfig {
                slab_size: config.slab_size,
                validate: true,
                ..FastqConfig::default()
            },
        );
        consume_fastq(&mut reader)
    })
}

#[cfg(all(feature = "bgzf", feature = "libdeflate"))]
fn measure_bgzf_libdeflate_serial_trusted_pack(
    name: &str,
    input: &[u8],
    config: &Config,
) -> Result<Measurement> {
    measure(name, input.len(), config.iters, || {
        let source = microraptor::BgzfReader::with_inflate_backend(
            input,
            microraptor::BgzfInflateBackend::Libdeflate,
        );
        consume_trusted_fastq_read_with_pack(source, fastq_config(config))
    })
}

#[cfg(all(feature = "bgzf", feature = "libdeflate"))]
fn measure_bgzf_libdeflate_adaptive_trusted_pack(
    name: &str,
    input: &[u8],
    config: &Config,
) -> Result<Measurement> {
    let owned: Arc<[u8]> = Arc::from(input);
    measure_with_bgzf_metrics(name, input.len(), config.iters, |metrics| {
        let source = microraptor::BgzfAutoReader::with_config(
            std::io::Cursor::new(Arc::clone(&owned)),
            input.len() as u64,
            bgzf_config(config)
                .with_inflate_backend(microraptor::BgzfInflateBackend::Libdeflate)
                .with_metrics(metrics),
        )?;
        consume_trusted_fastq_read_with_pack(source, fastq_config(config))
    })
}

fn measure_pack(name: &str, input: &[u8], config: &Config) -> Result<Measurement> {
    measure(name, input.len(), config.iters, || {
        let mut reader = FastqReader::with_config(
            std::io::Cursor::new(input),
            FastqConfig {
                slab_size: config.slab_size,
                validate: true,
                ..FastqConfig::default()
            },
        );
        consume_fastq_with_pack(&mut reader)
    })
}

fn measure_trusted_pack(name: &str, input: &[u8], config: &Config) -> Result<Measurement> {
    measure(name, input.len(), config.iters, || {
        consume_trusted_fastq_read_with_pack(
            std::io::Cursor::new(input),
            FastqConfig {
                slab_size: config.slab_size,
                validate: true,
                ..FastqConfig::default()
            },
        )
    })
}

fn measure_direct_pack(name: &str, input: &[u8], config: &Config) -> Result<Measurement> {
    measure(name, input.len(), config.iters, || {
        consume_trusted_fastq_read_direct_with_pack(
            std::io::Cursor::new(input),
            FastqConfig {
                slab_size: config.slab_size,
                validate: true,
                ..FastqConfig::default()
            },
        )
    })
}

fn measure_path_fastq(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        let mut reader =
            FastqReader::with_config(open_bench_read(path, config)?, fastq_config(config));
        consume_fastq(&mut reader)
    })
}

#[cfg(feature = "mmap")]
fn measure_path_fastq_mmap(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        let mut stats = StreamStats::default();
        microraptor::visit_fastq_mmap(path, fastq_config(config), |record| {
            stats.observe_record(record.name(), record.seq(), record.qual());
            Ok(())
        })?;
        Ok(stats)
    })
}

fn measure_path_fasta(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        let mut reader = FastaReader::with_config(
            open_bench_read(path, config)?,
            microraptor::FastaConfig {
                batch_records: fasta_batch_records(config),
                ..microraptor::FastaConfig::default()
            },
        );
        consume_fasta(&mut reader)
    })
}

#[cfg(feature = "mmap")]
fn measure_path_fasta_mmap(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        let mut stats = StreamStats::default();
        microraptor::visit_fasta_mmap(path, |record| {
            stats.observe_sequence_record(record.name(), record.seq());
            Ok(())
        })?;
        Ok(stats)
    })
}

fn measure_path_pack(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        let mut reader =
            FastqReader::with_config(open_bench_read(path, config)?, fastq_config(config));
        consume_fastq_with_pack(&mut reader)
    })
}

fn measure_path_trusted_pack(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    #[cfg(feature = "bgzf")]
    if path_has_bgzf_header(path)? {
        return measure_path_bgzf_adaptive_trusted_pack(name, path, input_bytes, config);
    }

    measure(name, input_bytes, config.iters, || {
        consume_trusted_fastq_read_with_pack(open_bench_read(path, config)?, fastq_config(config))
    })
}

fn measure_path_direct_pack(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        consume_trusted_fastq_read_direct_with_pack(
            open_bench_read(path, config)?,
            fastq_config(config),
        )
    })
}

fn measure_paired_path_fastq(
    name: &str,
    first: &Path,
    second: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        let mut reader = microraptor::PairedFastqReader::from_fastq_readers(
            FastqReader::with_config(open_bench_read(first, config)?, fastq_config(config)),
            FastqReader::with_config(open_bench_read(second, config)?, fastq_config(config)),
        );
        consume_paired_fastq(&mut reader)
    })
}

fn measure_paired_path_pack(
    name: &str,
    first: &Path,
    second: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        let mut reader = microraptor::PairedFastqReader::from_fastq_readers(
            FastqReader::with_config(open_bench_read(first, config)?, fastq_config(config)),
            FastqReader::with_config(open_bench_read(second, config)?, fastq_config(config)),
        );
        let mut ctx = PackContext::default();
        while let Some(batch) = reader.next_pair_batch()? {
            for pair in batch.pairs() {
                let first = pair.first();
                ctx.observe_packed(first.name(), first.seq(), first.qual())?;
                let second = pair.second();
                ctx.observe_packed(second.name(), second.seq(), second.qual())?;
            }
        }
        Ok(ctx.stats)
    })
}

fn measure_paired_path_trusted_pack(
    name: &str,
    first: &Path,
    second: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        let mut stats = StreamStats::default();
        microraptor::pack::pack_trusted_paired_fastq_read(
            open_bench_read(first, config)?,
            open_bench_read(second, config)?,
            fastq_config(config),
            PairValidation::FastSlash,
            |pair| {
                observe_trusted_record(&mut stats, pair.first);
                observe_trusted_record(&mut stats, pair.second);
                Ok(())
            },
        )?;
        Ok(stats)
    })
}

fn observe_trusted_record(stats: &mut StreamStats, record: TrustedPackedRecord<'_>) {
    stats.observe_record(record.name, record.seq, record.qual);
    stats.checksum = stats
        .checksum
        .wrapping_add(record.summary.bases.canonical_bases() as u64)
        .wrapping_add(record.summary.qualities.sum_phred);
}

#[derive(Default)]
struct PackContext {
    stats: StreamStats,
    packed: Vec<u8>,
    mask: Vec<u8>,
}

impl PackContext {
    fn observe_packed(&mut self, name: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
        let summary =
            pack_bases_and_summarize_qualities_into(seq, qual, &mut self.packed, &mut self.mask)
                .map_err(|e| microraptor::FastqError::Format(e.to_string()))?;
        self.stats.observe_record(name, seq, qual);
        self.stats.checksum = self
            .stats
            .checksum
            .wrapping_add(summary.bases.canonical_bases() as u64)
            .wrapping_add(summary.qualities.sum_phred);
        Ok(())
    }
}

fn consume_fastq_with_pack<R: std::io::Read>(reader: &mut FastqReader<R>) -> Result<StreamStats> {
    let mut ctx = PackContext::default();
    while let Some(batch) = reader.next_batch()? {
        for record in batch.records() {
            ctx.observe_packed(record.name(), record.seq(), record.qual())?;
        }
    }
    Ok(ctx.stats)
}

fn consume_paired_fastq<R1: std::io::Read, R2: std::io::Read>(
    reader: &mut microraptor::PairedFastqReader<R1, R2>,
) -> Result<StreamStats> {
    let mut stats = StreamStats::default();
    while let Some(batch) = reader.next_pair_batch()? {
        for pair in batch.pairs() {
            let first = pair.first();
            stats.observe_record(first.name(), first.seq(), first.qual());
            let second = pair.second();
            stats.observe_record(second.name(), second.seq(), second.qual());
        }
    }
    Ok(stats)
}

#[cfg(all(feature = "bgzf", feature = "libdeflate"))]
fn measure_path_bgzf_libdeflate_serial(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        let file = std::fs::File::open(path)?;
        let source = microraptor::BgzfReader::with_inflate_backend(
            file,
            microraptor::BgzfInflateBackend::Libdeflate,
        );
        let mut reader = FastqReader::with_config(source, fastq_config(config));
        consume_fastq(&mut reader)
    })
}

#[cfg(all(feature = "bgzf", feature = "libdeflate"))]
fn measure_path_bgzf_libdeflate_parallel(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        let file = std::fs::File::open(path)?;
        let source = microraptor::BgzfParallelReader::with_inflate_backend(
            file,
            config.workers,
            microraptor::BgzfInflateBackend::Libdeflate,
        )?;
        let mut reader = FastqReader::with_config(source, fastq_config(config));
        consume_fastq(&mut reader)
    })
}

#[cfg(feature = "bgzf")]
fn measure_path_bgzf_adaptive_trusted_pack(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure_with_bgzf_metrics(name, input_bytes, config.iters, |metrics| {
        let file = std::fs::File::open(path)?;
        let source = microraptor::BgzfAutoReader::with_config(
            file,
            input_bytes as u64,
            bgzf_config(config).with_metrics(metrics),
        )?;
        consume_trusted_fastq_read_with_pack(source, fastq_config(config))
    })
}

#[cfg(feature = "bgzf")]
fn measure_path_bgzf_direct_parallel_trusted_pack(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure_with_bgzf_metrics(name, input_bytes, config.iters, |metrics| {
        let file = std::fs::File::open(path)?;
        let source = microraptor::BgzfParallelReader::with_config(
            file,
            bgzf_config(config).with_metrics(metrics),
        )?;
        consume_trusted_fastq_read_with_pack(source, fastq_config(config))
    })
}

#[cfg(all(feature = "bgzf", feature = "libdeflate"))]
fn measure_path_bgzf_libdeflate_serial_trusted_pack(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        let file = std::fs::File::open(path)?;
        let source = microraptor::BgzfReader::with_inflate_backend(
            file,
            microraptor::BgzfInflateBackend::Libdeflate,
        );
        consume_trusted_fastq_read_with_pack(source, fastq_config(config))
    })
}

#[cfg(all(feature = "bgzf", feature = "libdeflate"))]
fn measure_path_bgzf_libdeflate_adaptive_trusted_pack(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure_with_bgzf_metrics(name, input_bytes, config.iters, |metrics| {
        let file = std::fs::File::open(path)?;
        let source = microraptor::BgzfAutoReader::with_config(
            file,
            input_bytes as u64,
            bgzf_config(config)
                .with_inflate_backend(microraptor::BgzfInflateBackend::Libdeflate)
                .with_metrics(metrics),
        )?;
        consume_trusted_fastq_read_with_pack(source, fastq_config(config))
    })
}

#[cfg(feature = "bgzf")]
fn path_has_bgzf_header(path: &Path) -> Result<bool> {
    Ok(matches!(
        microraptor::detect_file_input_kind(path)?,
        DetectedInputKind::Bgzf
    ))
}

fn fastq_config(config: &Config) -> FastqConfig {
    FastqConfig {
        slab_size: config.slab_size,
        validate: true,
        pair_validation: PairValidation::FastSlash,
        ..FastqConfig::default()
    }
}

fn fasta_batch_records(config: &Config) -> usize {
    (config.slab_size / 256).max(1)
}

#[cfg(feature = "bgzf")]
fn bgzf_config(config: &Config) -> microraptor::BgzfParallelConfig {
    microraptor::BgzfParallelConfig::new(config.workers)
        .with_parallel_min_compressed_bytes(config.bgzf_parallel_min_bytes)
}

#[cfg(feature = "bgzf")]
fn measure_with_bgzf_metrics<F>(
    name: &str,
    bytes: usize,
    iters: usize,
    mut f: F,
) -> Result<Measurement>
where
    F: FnMut(Arc<microraptor::BgzfPipelineMetrics>) -> Result<StreamStats>,
{
    let metrics = Arc::new(microraptor::BgzfPipelineMetrics::default());
    let mut measurement = measure(name, bytes, iters, || f(Arc::clone(&metrics)))?;
    let snapshot = metrics.snapshot();
    measurement
        .extras
        .push(("bgzf_job_queue_full", snapshot.job_queue_full));
    measurement
        .extras
        .push(("bgzf_result_queue_full", snapshot.result_queue_full));
    Ok(measurement)
}

fn measure<F>(name: &str, bytes: usize, iters: usize, mut f: F) -> Result<Measurement>
where
    F: FnMut() -> Result<StreamStats>,
{
    let mut best = Duration::MAX;
    let mut samples = Vec::with_capacity(iters.max(1));
    let mut last = StreamStats::default();
    for _ in 0..iters.max(1) {
        let start = Instant::now();
        last = f()?;
        let elapsed = start.elapsed();
        best = best.min(elapsed);
        samples.push(elapsed);
    }
    Ok(Measurement {
        name: name.to_string(),
        bytes,
        records: last.records,
        bases: last.bases,
        best,
        samples,
        checksum: last.checksum,
        extras: Vec::new(),
    })
}

#[cfg(feature = "gzip")]
fn gzip_bytes(raw: &[u8]) -> Result<Vec<u8>> {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    encoder.write_all(raw)?;
    Ok(encoder.finish()?)
}

fn emit_report(config: &Config, source: &str, input_bytes: usize, rows: &[Measurement]) {
    let report_rows = report_rows(rows);
    if config.json {
        println!(
            "{}",
            report::render_json(&report_config(config), source, input_bytes, &report_rows)
        );
    } else {
        report::print_table(source, input_bytes, &report_rows);
    }
}

fn report_config(config: &Config) -> report::ReportConfig<'_> {
    report::ReportConfig {
        mode: config.mode.as_str(),
        format: config.format.as_str(),
        records: config.records,
        read_len: config.read_len,
        iters: config.iters,
        slab_size: config.slab_size,
        workers: config.workers,
    }
}

fn report_rows(rows: &[Measurement]) -> Vec<report::ReportRow<'_>> {
    rows.iter()
        .map(|row| report::ReportRow {
            name: &row.name,
            bytes: row.bytes,
            records: row.records,
            bases: row.bases,
            best: row.best,
            samples: &row.samples,
            checksum: row.checksum,
            extras: &row.extras,
        })
        .collect()
}

fn parse_args() -> Config {
    let mut config = Config::default();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--records" => config.records = parse_next(&mut args, "--records"),
            "--read-len" => config.read_len = parse_next(&mut args, "--read-len"),
            "--iters" => config.iters = parse_next(&mut args, "--iters"),
            "--slab-size" => config.slab_size = parse_next(&mut args, "--slab-size"),
            "--workers" => config.workers = parse_next(&mut args, "--workers"),
            "--bgzf-parallel-min-bytes" => {
                config.bgzf_parallel_min_bytes = parse_next(&mut args, "--bgzf-parallel-min-bytes")
            }
            "--input" => config.input = Some(parse_path(&mut args, "--input")),
            "--format" => config.format = parse_format(&mut args, "--format"),
            "--paired-inputs" => {
                let first = parse_path(&mut args, "--paired-inputs");
                let second = parse_path(&mut args, "--paired-inputs");
                config.paired_inputs = Some((first, second));
            }
            "--mode" => config.mode = parse_mode(&mut args, "--mode"),
            "--check-bgzf-pack-regression" => {
                config
                    .bgzf_pack_check
                    .get_or_insert_with(BgzfPackCheck::default);
            }
            "--check-label" => {
                config
                    .bgzf_pack_check
                    .get_or_insert_with(BgzfPackCheck::default)
                    .label = parse_string(&mut args, "--check-label");
            }
            "--min-input-bytes" => {
                config
                    .bgzf_pack_check
                    .get_or_insert_with(BgzfPackCheck::default)
                    .min_input_bytes = parse_next(&mut args, "--min-input-bytes");
            }
            "--tolerance-pct" => {
                config
                    .bgzf_pack_check
                    .get_or_insert_with(BgzfPackCheck::default)
                    .tolerance_pct = parse_next::<u128>(&mut args, "--tolerance-pct");
            }
            "--skip-timing-checks" => {
                config
                    .bgzf_pack_check
                    .get_or_insert_with(BgzfPackCheck::default)
                    .check_timing = false;
            }
            "--profile-bgzf-parallel" => config.profile_bgzf_parallel = true,
            "--mmap" => config.mmap = true,
            "--json" => config.json = true,
            "--help" | "-h" => {
                print_help();
                std::process::exit(0);
            }
            other => {
                eprintln!("unknown argument: {other}");
                print_help();
                std::process::exit(2);
            }
        }
    }
    if config.input.is_some() && config.paired_inputs.is_some() {
        eprintln!("--input and --paired-inputs are mutually exclusive");
        std::process::exit(2);
    }
    if config.format == InputFormat::Fasta {
        if config.paired_inputs.is_some() {
            eprintln!("--format fasta does not support --paired-inputs");
            std::process::exit(2);
        }
        if config.mode.includes_pack() {
            eprintln!("--format fasta supports --mode parse only");
            std::process::exit(2);
        }
    }
    config
}

fn parse_next<T>(args: &mut impl Iterator<Item = String>, flag: &str) -> T
where
    T: std::str::FromStr,
{
    let Some(value) = args.next() else {
        eprintln!("{flag} requires a value");
        std::process::exit(2);
    };
    match value.parse() {
        Ok(v) => v,
        Err(_) => {
            eprintln!("{flag} requires a numeric value, got {value}");
            std::process::exit(2);
        }
    }
}

fn parse_path(args: &mut impl Iterator<Item = String>, flag: &str) -> PathBuf {
    let Some(value) = args.next() else {
        eprintln!("{flag} requires a path");
        std::process::exit(2);
    };
    PathBuf::from(value)
}

fn parse_string(args: &mut impl Iterator<Item = String>, flag: &str) -> String {
    let Some(value) = args.next() else {
        eprintln!("{flag} requires a value");
        std::process::exit(2);
    };
    value
}

fn parse_mode(args: &mut impl Iterator<Item = String>, flag: &str) -> Mode {
    let Some(value) = args.next() else {
        eprintln!("{flag} requires one of: all, parse, pack");
        std::process::exit(2);
    };
    match value.as_str() {
        "all" => Mode::All,
        "parse" => Mode::Parse,
        "pack" => Mode::Pack,
        _ => {
            eprintln!("{flag} requires one of: all, parse, pack; got {value}");
            std::process::exit(2);
        }
    }
}

fn parse_format(args: &mut impl Iterator<Item = String>, flag: &str) -> InputFormat {
    let Some(value) = args.next() else {
        eprintln!("{flag} requires one of: fastq, fasta");
        std::process::exit(2);
    };
    match value.as_str() {
        "fastq" => InputFormat::Fastq,
        "fasta" => InputFormat::Fasta,
        _ => {
            eprintln!("{flag} requires one of: fastq, fasta; got {value}");
            std::process::exit(2);
        }
    }
}

fn print_help() {
    eprintln!(
        "microraptor-bench [--input PATH | --paired-inputs R1 R2] [--format fastq|fasta] [--mode all|parse|pack] [--records N] [--read-len N] [--iters N] [--slab-size BYTES] [--workers N] [--bgzf-parallel-min-bytes N] [--json] [--mmap] [--check-bgzf-pack-regression] [--check-label NAME] [--min-input-bytes N] [--tolerance-pct N] [--skip-timing-checks] [--profile-bgzf-parallel]"
    );
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    #[test]
    fn real_input_benchmark_reads_path() {
        let path = std::env::temp_dir().join(format!(
            "microraptor-bench-real-input-{}.fastq",
            std::process::id()
        ));
        let mut file = std::fs::File::create(&path).unwrap();
        file.write_all(b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nJJJJ\n")
            .unwrap();
        drop(file);

        let config = Config {
            iters: 1,
            input: Some(path.clone()),
            ..Config::default()
        };
        let result = run_real_input(&path, &config);

        std::fs::remove_file(path).unwrap();
        result.unwrap();
    }

    #[test]
    fn paired_input_benchmark_reads_paths() {
        let r1_path =
            std::env::temp_dir().join(format!("microraptor-bench-r1-{}.fastq", std::process::id()));
        let r2_path =
            std::env::temp_dir().join(format!("microraptor-bench-r2-{}.fastq", std::process::id()));
        std::fs::write(&r1_path, b"@frag/1\nACGT\n+\nIIII\n").unwrap();
        std::fs::write(&r2_path, b"@frag/2\nTGCA\n+\nJJJJ\n").unwrap();

        let config = Config {
            iters: 1,
            paired_inputs: Some((r1_path.clone(), r2_path.clone())),
            ..Config::default()
        };
        let result = run_paired_input(&r1_path, &r2_path, &config);

        std::fs::remove_file(r1_path).unwrap();
        std::fs::remove_file(r2_path).unwrap();
        result.unwrap();
    }
}
