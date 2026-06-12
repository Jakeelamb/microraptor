use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use microraptor::Result;

#[derive(Debug)]
struct Config {
    out_dir: PathBuf,
    records: usize,
    read_len: usize,
    pattern: Pattern,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            out_dir: PathBuf::from("target/bench-inputs"),
            records: 100_000,
            read_len: 150,
            pattern: Pattern::Cyclic,
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum Pattern {
    Cyclic,
    Entropy,
}

fn main() {
    if let Err(err) = run() {
        eprintln!("{err}");
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let config = parse_args();
    fs::create_dir_all(&config.out_dir)?;

    let single = build_single_end(config.records, config.read_len, config.pattern);
    let interleaved = build_interleaved(config.records, config.read_len, config.pattern);
    let (r1, r2) = build_paired(config.records, config.read_len, config.pattern);

    write_file(config.out_dir.join("single.fastq"), &single)?;
    write_file(config.out_dir.join("interleaved.fastq"), &interleaved)?;
    write_file(config.out_dir.join("r1.fastq"), &r1)?;
    write_file(config.out_dir.join("r2.fastq"), &r2)?;

    #[cfg(feature = "gzip")]
    {
        write_gzip(config.out_dir.join("single.fastq.gz"), &single)?;
        write_gzip(config.out_dir.join("interleaved.fastq.gz"), &interleaved)?;
        write_gzip(config.out_dir.join("r1.fastq.gz"), &r1)?;
        write_gzip(config.out_dir.join("r2.fastq.gz"), &r2)?;
    }

    #[cfg(feature = "bgzf")]
    {
        write_bgzf(config.out_dir.join("single.fastq.bgz"), &single)?;
        write_bgzf(config.out_dir.join("interleaved.fastq.bgz"), &interleaved)?;
        write_bgzf(config.out_dir.join("r1.fastq.bgz"), &r1)?;
        write_bgzf(config.out_dir.join("r2.fastq.bgz"), &r2)?;
    }

    println!("wrote fixtures to {}", config.out_dir.display());
    Ok(())
}

fn build_single_end(records: usize, read_len: usize, pattern: Pattern) -> Vec<u8> {
    let mut out = Vec::with_capacity(records.saturating_mul(read_len + 32));
    for i in 0..records {
        push_record(&mut out, b"r", i, None, read_len, 0, pattern);
    }
    out
}

fn build_interleaved(pairs: usize, read_len: usize, pattern: Pattern) -> Vec<u8> {
    let mut out = Vec::with_capacity(pairs.saturating_mul((read_len + 36) * 2));
    for i in 0..pairs {
        push_record(&mut out, b"frag", i, Some(1), read_len, 0, pattern);
        push_record(&mut out, b"frag", i, Some(2), read_len, 1, pattern);
    }
    out
}

fn build_paired(pairs: usize, read_len: usize, pattern: Pattern) -> (Vec<u8>, Vec<u8>) {
    let mut r1 = Vec::with_capacity(pairs.saturating_mul(read_len + 36));
    let mut r2 = Vec::with_capacity(pairs.saturating_mul(read_len + 36));
    for i in 0..pairs {
        push_record(&mut r1, b"frag", i, Some(1), read_len, 0, pattern);
        push_record(&mut r2, b"frag", i, Some(2), read_len, 1, pattern);
    }
    (r1, r2)
}

fn push_record(
    out: &mut Vec<u8>,
    prefix: &[u8],
    index: usize,
    mate: Option<u8>,
    read_len: usize,
    phase: usize,
    pattern: Pattern,
) {
    out.push(b'@');
    out.extend_from_slice(prefix);
    push_usize_decimal(index, out);
    if let Some(mate) = mate {
        out.push(b'/');
        out.push(b'0' + mate);
    }
    out.push(b'\n');

    push_bases(out, index, phase, read_len, pattern);
    out.extend_from_slice(b"\n+\n");
    push_qualities(out, index, phase, read_len, pattern);
    out.push(b'\n');
}

fn push_bases(out: &mut Vec<u8>, index: usize, phase: usize, read_len: usize, pattern: Pattern) {
    let bases = b"ACGT";
    match pattern {
        Pattern::Cyclic => {
            for j in 0..read_len {
                out.push(bases[(index + j + phase) & 3]);
            }
        }
        Pattern::Entropy => {
            let mut state = rng_seed(index, phase, 0xa076_1d64_78bd_642f);
            for _ in 0..read_len {
                out.push(bases[(next_u64(&mut state) as usize) & 3]);
            }
        }
    }
}

fn push_qualities(
    out: &mut Vec<u8>,
    index: usize,
    phase: usize,
    read_len: usize,
    pattern: Pattern,
) {
    match pattern {
        Pattern::Cyclic => out.extend(std::iter::repeat_n(b'I', read_len)),
        Pattern::Entropy => {
            let mut state = rng_seed(index, phase, 0xe703_7ed1_a0b4_28db);
            for _ in 0..read_len {
                out.push(33 + (next_u64(&mut state) % 41) as u8);
            }
        }
    }
}

fn rng_seed(index: usize, phase: usize, salt: u64) -> u64 {
    salt ^ ((index as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15)) ^ ((phase as u64) << 32)
}

fn next_u64(state: &mut u64) -> u64 {
    *state ^= *state >> 12;
    *state ^= *state << 25;
    *state ^= *state >> 27;
    state.wrapping_mul(0x2545_f491_4f6c_dd1d)
}

fn push_usize_decimal(mut n: usize, out: &mut Vec<u8>) {
    if n == 0 {
        out.push(b'0');
        return;
    }
    let mut buf = [0_u8; 20];
    let mut i = buf.len();
    while n != 0 {
        i -= 1;
        buf[i] = b'0' + (n % 10) as u8;
        n /= 10;
    }
    out.extend_from_slice(&buf[i..]);
}

fn write_file(path: PathBuf, bytes: &[u8]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writer.write_all(bytes)?;
    writer.flush()?;
    Ok(())
}

#[cfg(feature = "gzip")]
fn write_gzip(path: PathBuf, bytes: &[u8]) -> Result<()> {
    let file = File::create(path)?;
    let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::fast());
    encoder.write_all(bytes)?;
    encoder.finish()?;
    Ok(())
}

#[cfg(feature = "bgzf")]
fn write_bgzf(path: PathBuf, bytes: &[u8]) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = microraptor::BgzfWriter::new(BufWriter::new(file));
    writer.write_all(bytes)?;
    writer.finish()?;
    Ok(())
}

fn parse_args() -> Config {
    let mut config = Config::default();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--out-dir" => config.out_dir = parse_path(&mut args, "--out-dir"),
            "--records" => config.records = parse_usize(&mut args, "--records"),
            "--read-len" => config.read_len = parse_usize(&mut args, "--read-len"),
            "--pattern" => config.pattern = parse_pattern(&mut args, "--pattern"),
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

fn parse_path(args: &mut impl Iterator<Item = String>, flag: &str) -> PathBuf {
    let Some(value) = args.next() else {
        eprintln!("{flag} requires a path");
        std::process::exit(2);
    };
    PathBuf::from(value)
}

fn parse_usize(args: &mut impl Iterator<Item = String>, flag: &str) -> usize {
    let Some(value) = args.next() else {
        eprintln!("{flag} requires a value");
        std::process::exit(2);
    };
    match value.parse() {
        Ok(value) => value,
        Err(_) => {
            eprintln!("{flag} requires an unsigned integer, got {value}");
            std::process::exit(2);
        }
    }
}

fn parse_pattern(args: &mut impl Iterator<Item = String>, flag: &str) -> Pattern {
    let Some(value) = args.next() else {
        eprintln!("{flag} requires one of: cyclic, entropy");
        std::process::exit(2);
    };
    match value.as_str() {
        "cyclic" => Pattern::Cyclic,
        "entropy" => Pattern::Entropy,
        _ => {
            eprintln!("{flag} requires one of: cyclic, entropy; got {value}");
            std::process::exit(2);
        }
    }
}

fn print_help() {
    eprintln!(
        "microraptor-fixture [--out-dir PATH] [--records N] [--read-len N] [--pattern cyclic|entropy]"
    );
}
