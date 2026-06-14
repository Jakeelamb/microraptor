use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use microraptor::Result;
use microraptor::benchutil::{
    SyntheticAlphabet as Alphabet, SyntheticFastaLayout as FastaLayout,
    SyntheticPattern as Pattern, synthetic_fasta_with_options, synthetic_fastq_with_pattern,
    synthetic_interleaved_fastq_with_pattern, synthetic_paired_fastq_with_pattern,
};

#[derive(Debug)]
struct Config {
    out_dir: PathBuf,
    records: usize,
    read_len: usize,
    pattern: Pattern,
    format: Format,
    fasta_layout: FastaLayout,
    alphabet: Alphabet,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            out_dir: PathBuf::from("target/bench-inputs"),
            records: 100_000,
            read_len: 150,
            pattern: Pattern::Cyclic,
            format: Format::Fastq,
            fasta_layout: FastaLayout::TwoLine,
            alphabet: Alphabet::Dna,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Format {
    Fastq,
    Fasta,
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

    if config.format == Format::Fasta {
        let single = synthetic_fasta_with_options(
            config.records,
            config.read_len,
            config.pattern,
            config.fasta_layout,
            config.alphabet,
        );
        write_file(config.out_dir.join("single.fasta"), &single)?;

        #[cfg(feature = "gzip")]
        {
            write_gzip(config.out_dir.join("single.fasta.gz"), &single)?;
        }

        #[cfg(feature = "bgzf")]
        {
            write_bgzf(config.out_dir.join("single.fasta.bgz"), &single)?;
        }

        println!("wrote fixtures to {}", config.out_dir.display());
        return Ok(());
    }

    let single = synthetic_fastq_with_pattern(config.records, config.read_len, config.pattern);
    let interleaved =
        synthetic_interleaved_fastq_with_pattern(config.records, config.read_len, config.pattern);
    let (r1, r2) =
        synthetic_paired_fastq_with_pattern(config.records, config.read_len, config.pattern);

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
            "--format" => config.format = parse_format(&mut args, "--format"),
            "--fasta-layout" => {
                config.fasta_layout = parse_fasta_layout(&mut args, "--fasta-layout")
            }
            "--alphabet" => config.alphabet = parse_alphabet(&mut args, "--alphabet"),
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

fn parse_format(args: &mut impl Iterator<Item = String>, flag: &str) -> Format {
    let Some(value) = args.next() else {
        eprintln!("{flag} requires one of: fastq, fasta");
        std::process::exit(2);
    };
    match value.as_str() {
        "fastq" => Format::Fastq,
        "fasta" => Format::Fasta,
        _ => {
            eprintln!("{flag} requires one of: fastq, fasta; got {value}");
            std::process::exit(2);
        }
    }
}

fn parse_fasta_layout(args: &mut impl Iterator<Item = String>, flag: &str) -> FastaLayout {
    let Some(value) = args.next() else {
        eprintln!("{flag} requires one of: two-line, wrapped:N");
        std::process::exit(2);
    };
    if value == "two-line" {
        return FastaLayout::TwoLine;
    }
    if let Some(width) = value.strip_prefix("wrapped:") {
        return match width.parse() {
            Ok(width) => FastaLayout::Wrapped { width },
            Err(_) => {
                eprintln!("{flag} wrapped width must be an unsigned integer; got {value}");
                std::process::exit(2);
            }
        };
    }
    eprintln!("{flag} requires one of: two-line, wrapped:N; got {value}");
    std::process::exit(2);
}

fn parse_alphabet(args: &mut impl Iterator<Item = String>, flag: &str) -> Alphabet {
    let Some(value) = args.next() else {
        eprintln!("{flag} requires one of: dna, protein");
        std::process::exit(2);
    };
    match value.as_str() {
        "dna" => Alphabet::Dna,
        "protein" => Alphabet::Protein,
        _ => {
            eprintln!("{flag} requires one of: dna, protein; got {value}");
            std::process::exit(2);
        }
    }
}

fn print_help() {
    eprintln!(
        "microraptor-fixture [--out-dir PATH] [--format fastq|fasta] [--records N] [--read-len N] [--pattern cyclic|entropy] [--fasta-layout two-line|wrapped:N] [--alphabet dna|protein]"
    );
}
