use std::fmt::Write as _;
#[cfg(feature = "gzip")]
use std::io::Write;
use std::path::{Path, PathBuf};
#[cfg(feature = "bgzf")]
use std::sync::Arc;
use std::time::{Duration, Instant};

use microraptor::benchutil::{StreamStats, consume_fastq, synthetic_fastq};
use microraptor::pack::{pack_bases_into, summarize_qualities};
use microraptor::{FastqConfig, FastqReader, Result, open_fastq_with_config};

#[derive(Debug, Clone)]
struct Config {
    records: usize,
    read_len: usize,
    iters: usize,
    slab_size: usize,
    workers: usize,
    json: bool,
    input: Option<PathBuf>,
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
    if let Some(path) = config.input.as_deref() {
        return run_real_input(path, &config);
    }

    let raw = synthetic_fastq(config.records, config.read_len);
    #[cfg_attr(not(any(feature = "gzip", feature = "bgzf")), allow(unused_mut))]
    let mut measurements = vec![
        measure_fastq("raw", &raw, &config)?,
        measure_pack("pack-seq-qual", &raw, &config)?,
    ];

    #[cfg(feature = "gzip")]
    {
        let gzip = gzip_bytes(&raw)?;
        measurements.push(measure_gzip("gzip", &gzip, &config)?);
    }

    #[cfg(feature = "bgzf")]
    {
        let bgzf = microraptor::compress_bgzf_parallel(&raw, config.workers)?;
        measurements.push(measure_bgzf_serial("bgzf-serial", &bgzf, &config)?);
        measurements.push(measure_bgzf_parallel("bgzf-parallel", &bgzf, &config)?);
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

fn run_real_input(path: &Path, config: &Config) -> Result<()> {
    let input_bytes = usize::try_from(std::fs::metadata(path)?.len())
        .map_err(|_| microraptor::FastqError::Format("input file is too large".into()))?;
    let measurements = vec![
        measure_path_fastq("file-auto", path, input_bytes, config)?,
        measure_path_pack("file-pack-seq-qual", path, input_bytes, config)?,
    ];
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
        let mut stats = StreamStats::default();
        let mut packed = Vec::new();
        let mut mask = Vec::new();
        while let Some(batch) = reader.next_batch()? {
            for record in batch.records() {
                let seq = record.seq();
                let qual = record.qual();
                let summary = pack_bases_into(seq, &mut packed, &mut mask);
                let q = summarize_qualities(qual)
                    .map_err(|e| microraptor::FastqError::Format(e.to_string()))?;
                stats.observe_record(record.name(), seq, qual);
                stats.checksum = stats
                    .checksum
                    .wrapping_add(summary.canonical_bases() as u64)
                    .wrapping_add(q.sum_phred);
            }
        }
        Ok(stats)
    })
}

fn measure_path_fastq(
    name: &str,
    path: &Path,
    input_bytes: usize,
    config: &Config,
) -> Result<Measurement> {
    measure(name, input_bytes, config.iters, || {
        let mut reader = open_fastq_with_config(path, fastq_config(config))?;
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
        let mut reader = open_fastq_with_config(path, fastq_config(config))?;
        let mut stats = StreamStats::default();
        let mut packed = Vec::new();
        let mut mask = Vec::new();
        while let Some(batch) = reader.next_batch()? {
            for record in batch.records() {
                let seq = record.seq();
                let qual = record.qual();
                let summary = pack_bases_into(seq, &mut packed, &mut mask);
                let q = summarize_qualities(qual)
                    .map_err(|e| microraptor::FastqError::Format(e.to_string()))?;
                stats.observe_record(record.name(), seq, qual);
                stats.checksum = stats
                    .checksum
                    .wrapping_add(summary.canonical_bases() as u64)
                    .wrapping_add(q.sum_phred);
            }
        }
        Ok(stats)
    })
}

fn fastq_config(config: &Config) -> FastqConfig {
    FastqConfig {
        slab_size: config.slab_size,
        validate: true,
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
        "{{\"source\":{},\"records\":{},\"read_len\":{},\"iters\":{},\"slab_size\":{},\"workers\":{},\"input_bytes\":{},\"measurements\":[",
        JsonStr(source),
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

fn print_help() {
    eprintln!(
        "microraptor-bench [--input PATH] [--records N] [--read-len N] [--iters N] [--slab-size BYTES] [--workers N] [--json]"
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
    fn json_string_escapes_control_characters() {
        assert_eq!(JsonStr("a\"b\\c\n").to_string(), "\"a\\\"b\\\\c\\n\"");
    }
}
