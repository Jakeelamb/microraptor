use std::fmt::Write as _;
#[cfg(feature = "gzip")]
use std::io::Write;
#[cfg(feature = "bgzf")]
use std::sync::Arc;
use std::time::{Duration, Instant};

use microraptor::benchutil::{StreamStats, consume_fastq, synthetic_fastq};
use microraptor::pack::{pack_bases_into, summarize_qualities};
use microraptor::{FastqConfig, FastqReader, Result};

#[derive(Debug, Clone)]
struct Config {
    records: usize,
    read_len: usize,
    iters: usize,
    slab_size: usize,
    workers: usize,
    json: bool,
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
        println!("{}", render_json(&config, raw.len(), &measurements));
    } else {
        print_table(raw.len(), &measurements);
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

fn print_table(raw_bytes: usize, rows: &[Measurement]) {
    println!("raw_bytes\t{raw_bytes}");
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

fn render_json(config: &Config, raw_bytes: usize, rows: &[Measurement]) -> String {
    let mut out = String::new();
    let _ = write!(
        out,
        "{{\"records\":{},\"read_len\":{},\"iters\":{},\"slab_size\":{},\"workers\":{},\"raw_bytes\":{},\"measurements\":[",
        config.records, config.read_len, config.iters, config.slab_size, config.workers, raw_bytes
    );
    for (i, row) in rows.iter().enumerate() {
        if i != 0 {
            out.push(',');
        }
        let nanos = row.best.as_nanos();
        let _ = write!(
            out,
            "{{\"name\":\"{}\",\"input_bytes\":{},\"records\":{},\"bases\":{},\"best_ns\":{},\"checksum\":{}}}",
            row.name, row.bytes, row.records, row.bases, nanos, row.checksum
        );
    }
    out.push_str("]}");
    out
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

fn print_help() {
    eprintln!(
        "microraptor-bench [--records N] [--read-len N] [--iters N] [--slab-size BYTES] [--workers N] [--json]"
    );
}
