use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use microraptor::Result;

#[derive(Debug)]
struct Config {
    out_dir: PathBuf,
    records: usize,
    read_len: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            out_dir: PathBuf::from("target/bench-inputs"),
            records: 100_000,
            read_len: 150,
        }
    }
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

    let single = build_single_end(config.records, config.read_len);
    let interleaved = build_interleaved(config.records, config.read_len);
    let (r1, r2) = build_paired(config.records, config.read_len);

    write_file(config.out_dir.join("single.fastq"), &single)?;
    write_file(config.out_dir.join("interleaved.fastq"), &interleaved)?;
    write_file(config.out_dir.join("r1.fastq"), &r1)?;
    write_file(config.out_dir.join("r2.fastq"), &r2)?;

    #[cfg(feature = "gzip")]
    {
        write_gzip(config.out_dir.join("single.fastq.gz"), &single)?;
        write_gzip(config.out_dir.join("interleaved.fastq.gz"), &interleaved)?;
    }

    #[cfg(feature = "bgzf")]
    {
        write_bgzf(config.out_dir.join("single.fastq.bgz"), &single)?;
        write_bgzf(config.out_dir.join("interleaved.fastq.bgz"), &interleaved)?;
    }

    println!("wrote fixtures to {}", config.out_dir.display());
    Ok(())
}

fn build_single_end(records: usize, read_len: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(records.saturating_mul(read_len + 32));
    for i in 0..records {
        push_record(&mut out, b"r", i, None, read_len, 0);
    }
    out
}

fn build_interleaved(pairs: usize, read_len: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(pairs.saturating_mul((read_len + 36) * 2));
    for i in 0..pairs {
        push_record(&mut out, b"frag", i, Some(1), read_len, 0);
        push_record(&mut out, b"frag", i, Some(2), read_len, 1);
    }
    out
}

fn build_paired(pairs: usize, read_len: usize) -> (Vec<u8>, Vec<u8>) {
    let mut r1 = Vec::with_capacity(pairs.saturating_mul(read_len + 36));
    let mut r2 = Vec::with_capacity(pairs.saturating_mul(read_len + 36));
    for i in 0..pairs {
        push_record(&mut r1, b"frag", i, Some(1), read_len, 0);
        push_record(&mut r2, b"frag", i, Some(2), read_len, 1);
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
) {
    out.push(b'@');
    out.extend_from_slice(prefix);
    push_usize_decimal(index, out);
    if let Some(mate) = mate {
        out.push(b'/');
        out.push(b'0' + mate);
    }
    out.push(b'\n');

    let bases = b"ACGT";
    for j in 0..read_len {
        out.push(bases[(index + j + phase) & 3]);
    }
    out.extend_from_slice(b"\n+\n");
    out.extend(std::iter::repeat_n(b'I', read_len));
    out.push(b'\n');
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

fn print_help() {
    eprintln!("microraptor-fixture [--out-dir PATH] [--records N] [--read-len N]");
}
