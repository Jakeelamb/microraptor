use std::fmt::Write as _;
use std::io::Read;
#[cfg(feature = "gzip")]
use std::io::Write;
use std::io::{Seek, SeekFrom};
use std::path::{Path, PathBuf};
#[cfg(feature = "bgzf")]
use std::sync::Arc;
use std::time::{Duration, Instant};

use microraptor::benchutil::{
    StreamStats, consume_fastq, consume_trusted_fastq_read_with_pack, synthetic_fastq,
};
use microraptor::pack::pack_bases_and_summarize_qualities_into;
use microraptor::{FastqConfig, FastqReader, PairValidation, Result};

enum BenchRead {
    Raw(std::fs::File),
    #[cfg(feature = "gzip")]
    Gzip(flate2::read::MultiGzDecoder<std::fs::File>),
    #[cfg(feature = "bgzf")]
    Bgzf(microraptor::BgzfReader<std::fs::File>),
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
    json: bool,
    input: Option<PathBuf>,
    paired_inputs: Option<(PathBuf, PathBuf)>,
    mode: Mode,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            records: 200_000,
            read_len: 150,
            iters: 5,
            slab_size: 8 * 1024 * 1024,
            workers: std::thread::available_parallelism().map_or(1, usize::from),
            json: false,
            input: None,
            paired_inputs: None,
            mode: Mode::All,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    All,
    Parse,
    Pack,
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
    checksum: u64,
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

    let raw = synthetic_fastq(config.records, config.read_len);
    #[cfg_attr(not(any(feature = "gzip", feature = "bgzf")), allow(unused_mut))]
    let mut measurements = Vec::new();
    if config.mode.includes_parse() {
        measurements.push(measure_fastq("raw", &raw, &config)?);
    }
    if config.mode.includes_pack() {
        measurements.push(measure_pack("pack-seq-qual", &raw, &config)?);
        measurements.push(measure_trusted_pack(
            "trusted-pack-seq-qual",
            &raw,
            &config,
        )?);
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
    }

    if config.json {
        println!(
            "{}",
            render_json(&config, "synthetic", raw.len(), &measurements)
        );
    } else {
        print_table("synthetic", raw.len(), &measurements);
    }
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
    }

    let source = format!("{}+{}", first.display(), second.display());
    if config.json {
        println!(
            "{}",
            render_json(config, &source, input_bytes, &measurements)
        );
    } else {
        print_table(&source, input_bytes, &measurements);
    }
    Ok(())
}

fn run_real_input(path: &Path, config: &Config) -> Result<()> {
    let input_bytes = checked_file_len(path)?;
    let mut measurements = Vec::new();
    if config.mode.includes_parse() {
        measurements.push(measure_path_fastq("file-auto", path, input_bytes, config)?);
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
        measurements.push(measure_path_pack(
            "file-pack-seq-qual",
            path,
            input_bytes,
            config,
        )?);
        measurements.push(measure_path_trusted_pack(
            "file-trusted-pack-seq-qual",
            path,
            input_bytes,
            config,
        )?);
    }
    let source = path.to_string_lossy();

    if config.json {
        println!(
            "{}",
            render_json(config, &source, input_bytes, &measurements)
        );
    } else {
        print_table(&source, input_bytes, &measurements);
    }
    Ok(())
}

fn checked_file_len(path: &Path) -> Result<usize> {
    usize::try_from(std::fs::metadata(path)?.len())
        .map_err(|_| microraptor::FastqError::Format("input file is too large".into()))
}

fn open_bench_read(path: &Path) -> Result<BenchRead> {
    let mut file = std::fs::File::open(path)?;
    let mut prefix = [0_u8; 18];
    let _n = file.read(&mut prefix)?;
    file.seek(SeekFrom::Start(0))?;

    #[cfg(feature = "bgzf")]
    if is_bgzf_header(&prefix[.._n]) {
        return Ok(BenchRead::Bgzf(microraptor::BgzfReader::new(file)));
    }

    #[cfg(feature = "gzip")]
    if _n >= 2 && prefix[..2] == [0x1f, 0x8b] {
        return Ok(BenchRead::Gzip(flate2::read::MultiGzDecoder::new(file)));
    }

    Ok(BenchRead::Raw(file))
}

#[cfg(feature = "bgzf")]
fn is_bgzf_header(prefix: &[u8]) -> bool {
    prefix.len() >= 18
        && prefix[0] == 31
        && prefix[1] == 139
        && prefix[2] == 8
        && prefix[3] & 4 != 0
        && u16::from_le_bytes([prefix[10], prefix[11]]) >= 6
        && prefix[12] == b'B'
        && prefix[13] == b'C'
        && prefix[14] == 2
        && prefix[15] == 0
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

fn measure_path_fastq(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        let mut reader = FastqReader::with_config(open_bench_read(path)?, fastq_config(config));
        consume_fastq(&mut reader)
    })
}

fn measure_path_pack(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        let mut reader = FastqReader::with_config(open_bench_read(path)?, fastq_config(config));
        consume_fastq_with_pack(&mut reader)
    })
}

fn measure_path_trusted_pack(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        consume_trusted_fastq_read_with_pack(open_bench_read(path)?, fastq_config(config))
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
            FastqReader::with_config(open_bench_read(first)?, fastq_config(config)),
            FastqReader::with_config(open_bench_read(second)?, fastq_config(config)),
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
            FastqReader::with_config(open_bench_read(first)?, fastq_config(config)),
            FastqReader::with_config(open_bench_read(second)?, fastq_config(config)),
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

#[cfg(all(feature = "bgzf", feature = "libdeflate"))]
fn path_has_bgzf_header(path: &Path) -> Result<bool> {
    let mut file = std::fs::File::open(path)?;
    let mut prefix = [0_u8; 18];
    let n = file.read(&mut prefix)?;
    Ok(n >= 18
        && prefix[0] == 31
        && prefix[1] == 139
        && prefix[2] == 8
        && prefix[3] & 4 != 0
        && u16::from_le_bytes([prefix[10], prefix[11]]) >= 6
        && prefix[12] == b'B'
        && prefix[13] == b'C'
        && prefix[14] == 2
        && prefix[15] == 0)
}

fn fastq_config(config: &Config) -> FastqConfig {
    FastqConfig {
        slab_size: config.slab_size,
        validate: true,
        pair_validation: PairValidation::FastSlash,
        ..FastqConfig::default()
    }
}

fn measure<F>(name: &str, bytes: usize, iters: usize, mut f: F) -> Result<Measurement>
where
    F: FnMut() -> Result<StreamStats>,
{
    let mut best = Duration::MAX;
    let mut last = StreamStats::default();
    for _ in 0..iters.max(1) {
        let start = Instant::now();
        last = f()?;
        best = best.min(start.elapsed());
    }
    Ok(Measurement {
        name: name.to_string(),
        bytes,
        records: last.records,
        bases: last.bases,
        best,
        checksum: last.checksum,
    })
}

#[cfg(feature = "gzip")]
fn gzip_bytes(raw: &[u8]) -> Result<Vec<u8>> {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    encoder.write_all(raw)?;
    Ok(encoder.finish()?)
}

fn print_table(source: &str, input_bytes: usize, rows: &[Measurement]) {
    println!("source\t{source}");
    println!("input_bytes\t{input_bytes}");
    println!("name\tinput_mb\tbest_ms\tinput_mib_s\trecords_s\tbases_s\tchecksum");
    for row in rows {
        let secs = row.best.as_secs_f64();
        let mib_s = (row.bytes as f64 / 1_048_576.0) / secs;
        let records_s = row.records as f64 / secs;
        let bases_s = row.bases as f64 / secs;
        println!(
            "{}\t{:.3}\t{:.3}\t{:.3}\t{:.0}\t{:.0}\t{}",
            row.name,
            row.bytes as f64 / 1_048_576.0,
            row.best.as_secs_f64() * 1000.0,
            mib_s,
            records_s,
            bases_s,
            row.checksum
        );
    }
}

fn render_json(config: &Config, source: &str, input_bytes: usize, rows: &[Measurement]) -> String {
    let mut out = String::new();
    let _ = write!(
        out,
        "{{\"source\":{},\"mode\":{},\"records\":{},\"read_len\":{},\"iters\":{},\"slab_size\":{},\"workers\":{},\"input_bytes\":{},\"measurements\":[",
        JsonStr(source),
        JsonStr(config.mode.as_str()),
        config.records,
        config.read_len,
        config.iters,
        config.slab_size,
        config.workers,
        input_bytes
    );
    for (i, row) in rows.iter().enumerate() {
        if i != 0 {
            out.push(',');
        }
        let nanos = row.best.as_nanos();
        let _ = write!(
            out,
            "{{\"name\":{},\"input_bytes\":{},\"records\":{},\"bases\":{},\"best_ns\":{},\"checksum\":{}}}",
            JsonStr(&row.name),
            row.bytes,
            row.records,
            row.bases,
            nanos,
            row.checksum
        );
    }
    out.push_str("]}");
    out
}

struct JsonStr<'a>(&'a str);

impl std::fmt::Display for JsonStr<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_char('"')?;
        for ch in self.0.chars() {
            match ch {
                '"' => f.write_str("\\\"")?,
                '\\' => f.write_str("\\\\")?,
                '\n' => f.write_str("\\n")?,
                '\r' => f.write_str("\\r")?,
                '\t' => f.write_str("\\t")?,
                ch if ch.is_control() => write!(f, "\\u{:04x}", ch as u32)?,
                ch => f.write_char(ch)?,
            }
        }
        f.write_char('"')
    }
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
            "--input" => config.input = Some(parse_path(&mut args, "--input")),
            "--paired-inputs" => {
                let first = parse_path(&mut args, "--paired-inputs");
                let second = parse_path(&mut args, "--paired-inputs");
                config.paired_inputs = Some((first, second));
            }
            "--mode" => config.mode = parse_mode(&mut args, "--mode"),
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
    config
}

fn parse_next(args: &mut impl Iterator<Item = String>, flag: &str) -> usize {
    let Some(value) = args.next() else {
        eprintln!("{flag} requires a value");
        std::process::exit(2);
    };
    match value.parse() {
        Ok(v) => v,
        Err(_) => {
            eprintln!("{flag} requires an unsigned integer, got {value}");
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

fn print_help() {
    eprintln!(
        "microraptor-bench [--input PATH | --paired-inputs R1 R2] [--mode all|parse|pack] [--records N] [--read-len N] [--iters N] [--slab-size BYTES] [--workers N] [--json]"
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

    #[test]
    fn json_string_escapes_control_characters() {
        assert_eq!(JsonStr("a\"b\\c\n").to_string(), "\"a\\\"b\\\\c\\n\"");
    }

    #[test]
    fn render_json_includes_mode() {
        let config = Config {
            mode: Mode::Pack,
            ..Config::default()
        };
        let json = render_json(&config, "synthetic", 0, &[]);
        assert!(json.contains("\"mode\":\"pack\""));
    }
}
