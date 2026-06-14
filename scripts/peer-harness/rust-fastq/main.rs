use std::env;
use std::fs;
use std::hint::black_box;
use std::io::{BufRead, BufReader, Cursor, Write};
use std::path::PathBuf;
use std::sync::OnceLock;
use std::time::{Duration, Instant};

use bio::io::fastq as bio_fastq;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression as GzipCompression;
use microraptor::benchutil::synthetic_fastq;
use microraptor::{visit_fastq_bytes, FastqConfig, FastqReader};
use noodles_fastq as noodles_fastq;
use seq_io::fastq::Record as SeqIoRecord;

type AppResult<T> = Result<T, Box<dyn std::error::Error>>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Stats {
    records: u64,
    bases: u64,
    checksum: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Consumer {
    FullChecksum,
    LightAccounting,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum InputCompression {
    Raw,
    Gzip,
}

impl InputCompression {
    fn as_str(self) -> &'static str {
        match self {
            Self::Raw => "raw",
            Self::Gzip => "gzip",
        }
    }
}

impl Stats {
    fn new() -> Self {
        Self {
            records: 0,
            bases: 0,
            checksum: 0xcbf2_9ce4_8422_2325,
        }
    }

    fn observe(&mut self, seq: &[u8], qual: &[u8]) {
        self.records += 1;
        self.bases += seq.len() as u64;
        match consumer() {
            Consumer::FullChecksum => {
                self.checksum = mix_bytes(self.checksum, seq);
                self.checksum = mix_bytes(self.checksum, qual);
            }
            Consumer::LightAccounting => {
                self.checksum = mix_record_shape(self.checksum, seq, qual);
            }
        }
    }
}

fn consumer() -> Consumer {
    static CONSUMER: OnceLock<Consumer> = OnceLock::new();
    *CONSUMER.get_or_init(|| match env::var("MICRORAPTOR_RUST_PEER_CONSUMER") {
        Ok(value) if value == "light" => Consumer::LightAccounting,
        Ok(value) if value == "full" => Consumer::FullChecksum,
        Ok(value) if !value.is_empty() => {
            panic!("unsupported MICRORAPTOR_RUST_PEER_CONSUMER={value}; expected full or light")
        }
        _ => Consumer::FullChecksum,
    })
}

fn microraptor_slab_size() -> Option<usize> {
    static SLAB_SIZE: OnceLock<Option<usize>> = OnceLock::new();
    *SLAB_SIZE.get_or_init(|| match env::var("MICRORAPTOR_RUST_PEER_SLAB_SIZE") {
        Ok(value) if !value.is_empty() => Some(
            value
                .parse()
                .unwrap_or_else(|_| panic!("invalid MICRORAPTOR_RUST_PEER_SLAB_SIZE={value}")),
        ),
        _ => None,
    })
}

fn microraptor_config(validate: bool) -> FastqConfig {
    let mut config = FastqConfig {
        validate,
        ..FastqConfig::default()
    };
    if let Some(slab_size) = microraptor_slab_size() {
        config.slab_size = slab_size;
    }
    config
}

#[derive(Debug)]
struct Row {
    tool: &'static str,
    stats: Stats,
    best: Duration,
}

#[derive(Debug)]
struct Config {
    input: Option<PathBuf>,
    records: usize,
    read_len: usize,
    iters: usize,
    out: PathBuf,
    compression: InputCompression,
}

fn main() -> AppResult<()> {
    let config = parse_args()?;
    let (source, input) = match config.input.as_ref() {
        Some(path) => (path.display().to_string(), fs::read(path)?),
        None => {
            let raw = synthetic_fastq(config.records, config.read_len);
            let input = match config.compression {
                InputCompression::Raw => raw,
                InputCompression::Gzip => gzip_bytes(&raw)?,
            };
            (
                format!(
                    "synthetic-{}:{}x{}",
                    config.compression.as_str(),
                    config.records,
                    config.read_len
                ),
                input,
            )
        }
    };

    let mut rows = vec![
        measure(
            "microraptor-stream",
            &input,
            config.iters,
            config.compression,
            parse_microraptor,
        )?,
        measure(
            "microraptor-stream-visitor",
            &input,
            config.iters,
            config.compression,
            parse_microraptor_visitor,
        )?,
    ];
    if config.compression == InputCompression::Raw {
        rows.insert(
            1,
            measure(
                "microraptor-slice-visitor",
                &input,
                config.iters,
                config.compression,
                parse_microraptor_slice_visitor,
            )?,
        );
    }
    if env::var_os("MICRORAPTOR_RUST_PEER_DIAGNOSTICS").is_some() {
        if config.compression == InputCompression::Raw {
            rows.push(measure(
                "microraptor-slice-visitor-no-validate",
                &input,
                config.iters,
                config.compression,
                parse_microraptor_slice_visitor_no_validate,
            )?);
        }
        rows.push(measure(
            "microraptor-visitor-no-validate",
            &input,
            config.iters,
            config.compression,
            parse_microraptor_visitor_no_validate,
        )?);
        rows.push(measure(
            "microraptor-no-validate",
            &input,
            config.iters,
            config.compression,
            parse_microraptor_no_validate,
        )?);
        rows.push(measure(
            "microraptor-record-refs",
            &input,
            config.iters,
            config.compression,
            parse_microraptor_record_refs,
        )?);
    }
    rows.extend([
        measure("seq_io", &input, config.iters, config.compression, parse_seq_io)?,
        measure(
            "noodles-fastq",
            &input,
            config.iters,
            config.compression,
            parse_noodles_fastq,
        )?,
        measure("bio", &input, config.iters, config.compression, parse_bio)?,
    ]);

    let reference = rows[0].stats;
    for row in &rows[1..] {
        if row.stats != reference {
            return Err(format!(
                "{} stats mismatch: {:?} != {:?}",
                row.tool, row.stats, reference
            )
            .into());
        }
    }

    let mut out = String::from(
        "tool\trecords\tbases\tbest_ms\trecords_s\tbases_s\tchecksum\tinput_bytes\titers\tsource\tcompression\n",
    );
    for row in rows {
        let ns = row.best.as_nanos().max(1);
        let records_s = row.stats.records as u128 * 1_000_000_000 / ns;
        let bases_s = row.stats.bases as u128 * 1_000_000_000 / ns;
        out.push_str(&format!(
            "{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            row.tool,
            row.stats.records,
            row.stats.bases,
            row.best.as_secs_f64() * 1000.0,
            records_s,
            bases_s,
            row.stats.checksum,
            input.len(),
            config.iters,
            source,
            config.compression.as_str()
        ));
    }
    fs::write(config.out, out)?;
    Ok(())
}

fn parse_args() -> AppResult<Config> {
    let mut input = None;
    let mut records = 100_000;
    let mut read_len = 150;
    let mut iters = 5;
    let mut out = PathBuf::from("rust-library-peers.tsv");
    let mut compression = InputCompression::Raw;

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--input" => input = Some(PathBuf::from(required_value(&mut args, "--input")?)),
            "--records" => records = required_value(&mut args, "--records")?.parse()?,
            "--read-len" => read_len = required_value(&mut args, "--read-len")?.parse()?,
            "--iters" => iters = required_value(&mut args, "--iters")?.parse()?,
            "--out" => out = PathBuf::from(required_value(&mut args, "--out")?),
            "--compression" => {
                compression = match required_value(&mut args, "--compression")?.as_str() {
                    "raw" => InputCompression::Raw,
                    "gzip" => InputCompression::Gzip,
                    value => return Err(format!("unsupported --compression {value}").into()),
                }
            }
            "--help" | "-h" => {
                print_help();
                std::process::exit(0);
            }
            _ => return Err(format!("unknown argument: {arg}").into()),
        }
    }

    Ok(Config {
        input,
        records,
        read_len,
        iters,
        out,
        compression,
    })
}

fn required_value(args: &mut impl Iterator<Item = String>, flag: &str) -> AppResult<String> {
    args.next()
        .ok_or_else(|| format!("missing value for {flag}").into())
}

fn print_help() {
    println!(
        "microraptor-rust-peer-bench [--input PATH] [--records N] [--read-len N] [--iters N] [--out PATH] [--compression raw|gzip]"
    );
}

fn measure(
    tool: &'static str,
    input: &[u8],
    iters: usize,
    compression: InputCompression,
    f: fn(&[u8], InputCompression) -> AppResult<Stats>,
) -> AppResult<Row> {
    let mut best = Duration::MAX;
    let mut stats = None;
    for _ in 0..iters {
        let start = Instant::now();
        let run_stats = f(input, compression)?;
        let elapsed = start.elapsed();
        black_box(run_stats.checksum);
        if let Some(previous) = stats {
            if previous != run_stats {
                return Err(format!("{tool} emitted unstable stats").into());
            }
        } else {
            stats = Some(run_stats);
        }
        best = best.min(elapsed);
    }

    Ok(Row {
        tool,
        stats: stats.expect("at least one iteration"),
        best,
    })
}

fn input_reader(input: &[u8], compression: InputCompression) -> Box<dyn BufRead + '_> {
    match compression {
        InputCompression::Raw => Box::new(BufReader::new(Cursor::new(input))),
        InputCompression::Gzip => Box::new(BufReader::new(MultiGzDecoder::new(Cursor::new(input)))),
    }
}

fn gzip_bytes(input: &[u8]) -> AppResult<Vec<u8>> {
    let mut encoder = GzEncoder::new(Vec::new(), GzipCompression::default());
    encoder.write_all(input)?;
    Ok(encoder.finish()?)
}

fn parse_microraptor(input: &[u8], compression: InputCompression) -> AppResult<Stats> {
    let mut reader = FastqReader::with_config(input_reader(input, compression), microraptor_config(true));
    let mut stats = Stats::new();
    while let Some(batch) = reader.next_batch()? {
        for record in batch.records() {
            stats.observe(record.seq(), record.qual());
        }
    }
    Ok(stats)
}

fn parse_microraptor_slice_visitor(input: &[u8], compression: InputCompression) -> AppResult<Stats> {
    if compression != InputCompression::Raw {
        return Err("slice visitor only supports raw resident FASTQ bytes".into());
    }
    let mut stats = Stats::new();
    visit_fastq_bytes(input, FastqConfig::default(), |record| {
        stats.observe(record.seq(), record.qual());
        Ok(())
    })?;
    Ok(stats)
}

fn parse_microraptor_slice_visitor_no_validate(input: &[u8], compression: InputCompression) -> AppResult<Stats> {
    if compression != InputCompression::Raw {
        return Err("slice visitor only supports raw resident FASTQ bytes".into());
    }
    let mut stats = Stats::new();
    visit_fastq_bytes(
        input,
        FastqConfig {
            validate: false,
            ..FastqConfig::default()
        },
        |record| {
            stats.observe(record.seq(), record.qual());
            Ok(())
        },
    )?;
    Ok(stats)
}

fn parse_microraptor_visitor(input: &[u8], compression: InputCompression) -> AppResult<Stats> {
    let mut reader = FastqReader::with_config(input_reader(input, compression), microraptor_config(true));
    let mut stats = Stats::new();
    reader.visit_records(|record| {
        stats.observe(record.seq(), record.qual());
        Ok(())
    })?;
    Ok(stats)
}

fn parse_microraptor_visitor_no_validate(input: &[u8], compression: InputCompression) -> AppResult<Stats> {
    let mut reader = FastqReader::with_config(input_reader(input, compression), microraptor_config(false));
    let mut stats = Stats::new();
    reader.visit_records(|record| {
        stats.observe(record.seq(), record.qual());
        Ok(())
    })?;
    Ok(stats)
}

fn parse_microraptor_no_validate(input: &[u8], compression: InputCompression) -> AppResult<Stats> {
    let mut reader = FastqReader::with_config(input_reader(input, compression), microraptor_config(false));
    let mut stats = Stats::new();
    while let Some(batch) = reader.next_batch()? {
        for record in batch.records() {
            stats.observe(record.seq(), record.qual());
        }
    }
    Ok(stats)
}

fn parse_microraptor_record_refs(input: &[u8], compression: InputCompression) -> AppResult<Stats> {
    let mut reader = FastqReader::with_config(input_reader(input, compression), microraptor_config(false));
    let mut stats = Stats::new();
    while let Some(batch) = reader.next_batch()? {
        let bytes = batch.bytes();
        for record in batch.record_refs() {
            let seq = &bytes[record.seq.start as usize..record.seq.end as usize];
            let qual = &bytes[record.qual.start as usize..record.qual.end as usize];
            stats.observe(seq, qual);
        }
    }
    Ok(stats)
}

fn parse_seq_io(input: &[u8], compression: InputCompression) -> AppResult<Stats> {
    let mut reader = seq_io::fastq::Reader::new(input_reader(input, compression));
    let mut stats = Stats::new();
    while let Some(record) = reader.next() {
        let record = record?;
        stats.observe(record.seq(), record.qual());
    }
    Ok(stats)
}

fn parse_noodles_fastq(input: &[u8], compression: InputCompression) -> AppResult<Stats> {
    let mut reader = noodles_fastq::io::Reader::new(input_reader(input, compression));
    let mut record = noodles_fastq::Record::default();
    let mut stats = Stats::new();
    loop {
        let n = reader.read_record(&mut record)?;
        if n == 0 {
            break;
        }
        stats.observe(record.sequence(), record.quality_scores());
    }
    Ok(stats)
}

fn parse_bio(input: &[u8], compression: InputCompression) -> AppResult<Stats> {
    let reader = bio_fastq::Reader::new(input_reader(input, compression));
    let mut stats = Stats::new();
    for record in reader.records() {
        let record = record?;
        stats.observe(record.seq(), record.qual());
    }
    Ok(stats)
}

fn mix_bytes(mut state: u64, bytes: &[u8]) -> u64 {
    for &byte in bytes {
        state ^= byte as u64;
        state = state.wrapping_mul(0x0000_0100_0000_01b3);
    }
    state
}

fn mix_record_shape(mut state: u64, seq: &[u8], qual: &[u8]) -> u64 {
    state ^= seq.len() as u64;
    state = state.rotate_left(5).wrapping_mul(0x0000_0100_0000_01b3);
    state ^= qual.len() as u64;
    state = state.rotate_left(7).wrapping_mul(0x0000_0100_0000_01b3);
    state ^= seq.first().copied().unwrap_or_default() as u64;
    state ^= (seq.last().copied().unwrap_or_default() as u64) << 8;
    state ^= (qual.first().copied().unwrap_or_default() as u64) << 16;
    state ^= (qual.last().copied().unwrap_or_default() as u64) << 24;
    state
}
