use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::ops::Range;
use std::path::PathBuf;

use microraptor::benchutil::{consume_fasta, consume_fastq};
use microraptor::{
    FastaIndex, FastaPartitionConfig, IndexedFastaReader, Result, build_fasta_index, open_fasta,
    open_fastq, plan_fasta_partitions,
};

#[cfg(feature = "bgzf")]
use microraptor::{
    BgzfIndexedFastaReader, DetectedInputKind, build_bgzf_index_strict, build_fasta_index_bgzf,
    detect_file_input_kind,
};

fn main() {
    if let Err(err) = run() {
        eprintln!("{err}");
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let mut args = std::env::args().skip(1);
    let Some(command) = args.next() else {
        print_help();
        return Ok(());
    };
    match command.as_str() {
        "stats" => stats(args.collect()),
        "checksum" => checksum(args.collect()),
        "fasta-index" => fasta_index(args.collect()),
        "fasta-fetch" => fasta_fetch(args.collect()),
        "fasta-partitions" => fasta_partitions(args.collect()),
        "fasta-chunks" => fasta_chunks(args.collect()),
        "verify-bgzf" => verify_bgzf(args.collect()),
        "--help" | "-h" | "help" => {
            print_help();
            Ok(())
        }
        _ => Err(microraptor::FastqError::Format(format!(
            "unknown command: {command}"
        ))),
    }
}

fn stats(args: Vec<String>) -> Result<()> {
    let mut format = String::from("fastq");
    let mut path = None;
    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--format" => {
                format = iter.next().ok_or_else(|| {
                    microraptor::FastqError::Format("--format requires a value".into())
                })?;
            }
            "--help" | "-h" => {
                print_help();
                return Ok(());
            }
            _ if path.is_none() => path = Some(PathBuf::from(arg)),
            _ => {
                return Err(microraptor::FastqError::Format(format!(
                    "unexpected stats argument: {arg}"
                )));
            }
        }
    }
    let path = required_path(path, "stats")?;
    match format.as_str() {
        "fastq" => {
            let mut reader = open_fastq(&path)?;
            let stats = consume_fastq(&mut reader)?;
            print_stream_stats(&stats);
        }
        "fasta" => {
            let mut reader = open_fasta(&path)?;
            let stats = consume_fasta(&mut reader)?;
            print_stream_stats(&stats);
        }
        _ => {
            return Err(microraptor::FastqError::Format(format!(
                "unsupported stats format: {format}"
            )));
        }
    }
    Ok(())
}

fn checksum(args: Vec<String>) -> Result<()> {
    let mut format = None;
    let mut path = None;
    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--format" => format = Some(required_value(&mut iter, "--format")?),
            "--help" | "-h" => {
                print_help();
                return Ok(());
            }
            _ if path.is_none() => path = Some(arg),
            _ => {
                return Err(microraptor::FastqError::Format(format!(
                    "unexpected checksum argument: {arg}"
                )));
            }
        }
    }
    let format = format
        .ok_or_else(|| microraptor::FastqError::Format("checksum requires --format".into()))?;
    let reader = input_reader(path.as_deref())?;
    let stats = match format.as_str() {
        "fastq" => {
            let mut reader = microraptor::FastqReader::new(reader);
            consume_fastq(&mut reader)?
        }
        "fasta" => {
            let mut reader = microraptor::FastaReader::new(reader);
            consume_fasta(&mut reader)?
        }
        "sam" => consume_sam_sequences(reader)?,
        _ => {
            return Err(microraptor::FastqError::Format(format!(
                "unsupported checksum format: {format}"
            )));
        }
    };
    print_stream_stats(&stats);
    Ok(())
}

fn input_reader(path: Option<&str>) -> Result<Box<dyn Read>> {
    match path {
        None | Some("-") => Ok(Box::new(std::io::stdin())),
        Some(path) => Ok(Box::new(File::open(path)?)),
    }
}

fn consume_sam_sequences<R: Read>(reader: R) -> Result<microraptor::benchutil::StreamStats> {
    let mut stats = microraptor::benchutil::StreamStats::default();
    let mut line = Vec::new();
    let mut reader = BufReader::new(reader);
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        let line = trim_line(&line);
        if line.is_empty() || line.starts_with(b"@") {
            continue;
        }
        let mut fields = line.split(|&b| b == b'\t');
        let name = fields
            .next()
            .ok_or_else(|| microraptor::FastqError::Format("SAM row missing QNAME".into()))?;
        let seq = fields
            .nth(8)
            .ok_or_else(|| microraptor::FastqError::Format("SAM row missing SEQ column".into()))?;
        if seq != b"*" {
            stats.observe_record(name, seq, b"");
        }
    }
    Ok(stats)
}

fn fasta_index(args: Vec<String>) -> Result<()> {
    let path = one_path(args, "fasta-index")?;
    #[cfg(feature = "bgzf")]
    let index = if detect_file_input_kind(&path)? == DetectedInputKind::Bgzf {
        build_fasta_index_bgzf(File::open(&path)?)?
    } else {
        build_fasta_index(File::open(&path)?)?
    };
    #[cfg(not(feature = "bgzf"))]
    let index = build_fasta_index(File::open(&path)?)?;
    print!("{}", index.to_fai_string());
    Ok(())
}

fn fasta_fetch(args: Vec<String>) -> Result<()> {
    let mut path = None;
    let mut fai = None;
    let mut name = None;
    let mut start = None;
    let mut end = None;
    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--fai" => fai = Some(PathBuf::from(required_value(&mut iter, "--fai")?)),
            "--name" => name = Some(required_value(&mut iter, "--name")?),
            "--start" => start = Some(parse_u64_value(&mut iter, "--start")?),
            "--end" => end = Some(parse_u64_value(&mut iter, "--end")?),
            "--help" | "-h" => {
                print_help();
                return Ok(());
            }
            _ if path.is_none() => path = Some(PathBuf::from(arg)),
            _ => {
                return Err(microraptor::FastqError::Format(format!(
                    "unexpected fasta-fetch argument: {arg}"
                )));
            }
        }
    }
    let path = required_path(path, "fasta-fetch")?;
    let fai = required_path(fai, "fasta-fetch --fai")?;
    let name =
        name.ok_or_else(|| microraptor::FastqError::Format("fasta-fetch requires --name".into()))?;
    let range = Range {
        start: start.ok_or_else(|| {
            microraptor::FastqError::Format("fasta-fetch requires --start".into())
        })?,
        end: end
            .ok_or_else(|| microraptor::FastqError::Format("fasta-fetch requires --end".into()))?,
    };
    let index = FastaIndex::from_fai_read(File::open(fai)?)?;

    #[cfg(feature = "bgzf")]
    let seq = if detect_file_input_kind(&path)? == DetectedInputKind::Bgzf {
        let bgzf_index = build_bgzf_index_strict(File::open(&path)?)?;
        let mut reader = BgzfIndexedFastaReader::new(File::open(&path)?, index, bgzf_index);
        reader.fetch(name.as_bytes(), range)?
    } else {
        let mut reader = IndexedFastaReader::new(File::open(&path)?, index);
        reader.fetch(name.as_bytes(), range)?
    };
    #[cfg(not(feature = "bgzf"))]
    let seq = {
        let mut reader = IndexedFastaReader::new(File::open(&path)?, index);
        reader.fetch(name.as_bytes(), range)?
    };

    std::io::stdout().write_all(&seq)?;
    std::io::stdout().write_all(b"\n")?;
    Ok(())
}

fn fasta_partitions(args: Vec<String>) -> Result<()> {
    let mut path = None;
    let mut fai = None;
    let mut parts = None;
    let mut overlap = 0_u64;
    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--fai" => fai = Some(PathBuf::from(required_value(&mut iter, "--fai")?)),
            "--parts" => parts = Some(parse_usize_value(&mut iter, "--parts")?),
            "--overlap" => overlap = parse_u64_value(&mut iter, "--overlap")?,
            "--help" | "-h" => {
                print_help();
                return Ok(());
            }
            _ if path.is_none() => path = Some(PathBuf::from(arg)),
            _ => {
                return Err(microraptor::FastqError::Format(format!(
                    "unexpected fasta-partitions argument: {arg}"
                )));
            }
        }
    }
    let _path = required_path(path, "fasta-partitions")?;
    let fai = required_path(fai, "fasta-partitions --fai")?;
    let parts = parts.ok_or_else(|| {
        microraptor::FastqError::Format("fasta-partitions requires --parts".into())
    })?;
    let index = FastaIndex::from_fai_read(File::open(fai)?)?;
    for partition in plan_fasta_partitions(&index, FastaPartitionConfig::new(parts, overlap))? {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            partition.partition_index,
            String::from_utf8_lossy(&partition.name),
            partition.core.start,
            partition.core.end,
            partition.fetch.start,
            partition.fetch.end,
            partition.core_offset_in_fetch()
        );
    }
    Ok(())
}

fn fasta_chunks(args: Vec<String>) -> Result<()> {
    let mut path = None;
    let mut fai = None;
    let mut name = None;
    let mut start = None;
    let mut end = None;
    let mut chunk_bases = None;
    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--fai" => fai = Some(PathBuf::from(required_value(&mut iter, "--fai")?)),
            "--name" => name = Some(required_value(&mut iter, "--name")?),
            "--start" => start = Some(parse_u64_value(&mut iter, "--start")?),
            "--end" => end = Some(parse_u64_value(&mut iter, "--end")?),
            "--chunk-bases" => chunk_bases = Some(parse_u64_value(&mut iter, "--chunk-bases")?),
            "--help" | "-h" => {
                print_help();
                return Ok(());
            }
            _ if path.is_none() => path = Some(PathBuf::from(arg)),
            _ => {
                return Err(microraptor::FastqError::Format(format!(
                    "unexpected fasta-chunks argument: {arg}"
                )));
            }
        }
    }
    let path = required_path(path, "fasta-chunks")?;
    let fai = required_path(fai, "fasta-chunks --fai")?;
    let name =
        name.ok_or_else(|| microraptor::FastqError::Format("fasta-chunks requires --name".into()))?;
    let range = Range {
        start: start.ok_or_else(|| {
            microraptor::FastqError::Format("fasta-chunks requires --start".into())
        })?,
        end: end
            .ok_or_else(|| microraptor::FastqError::Format("fasta-chunks requires --end".into()))?,
    };
    let chunk_bases = chunk_bases.ok_or_else(|| {
        microraptor::FastqError::Format("fasta-chunks requires --chunk-bases".into())
    })?;
    let index = FastaIndex::from_fai_read(File::open(fai)?)?;

    #[cfg(feature = "bgzf")]
    {
        if detect_file_input_kind(&path)? == DetectedInputKind::Bgzf {
            let bgzf_index = build_bgzf_index_strict(File::open(&path)?)?;
            let mut reader = BgzfIndexedFastaReader::new(File::open(&path)?, index, bgzf_index);
            for chunk in reader.reference_chunks(name.as_bytes(), range, chunk_bases)? {
                print_reference_chunk(&chunk?);
            }
            return Ok(());
        }
    }

    let mut reader = IndexedFastaReader::new(File::open(&path)?, index);
    for chunk in reader.reference_chunks(name.as_bytes(), range, chunk_bases)? {
        print_reference_chunk(&chunk?);
    }
    Ok(())
}

fn print_reference_chunk(chunk: &microraptor::FastaReferenceChunk) {
    println!(
        "{}\t{}\t{}",
        String::from_utf8_lossy(&chunk.name),
        chunk.global_offset,
        String::from_utf8_lossy(&chunk.seq)
    );
}

fn verify_bgzf(args: Vec<String>) -> Result<()> {
    let path = one_path(args, "verify-bgzf")?;
    #[cfg(feature = "bgzf")]
    {
        let index = build_bgzf_index_strict(File::open(path)?)?;
        println!("status\tok");
        println!("blocks\t{}", index.len());
        println!("uncompressed_len\t{}", index.uncompressed_len());
        println!("compressed_len\t{}", index.compressed_len());
        Ok(())
    }
    #[cfg(not(feature = "bgzf"))]
    {
        let _ = path;
        Err(microraptor::FastqError::Format(
            "verify-bgzf requires the bgzf feature".into(),
        ))
    }
}

fn print_stream_stats(stats: &microraptor::benchutil::StreamStats) {
    println!("records\t{}", stats.records);
    println!("bases\t{}", stats.bases);
    println!("qualities\t{}", stats.qualities);
    println!("name_bytes\t{}", stats.name_bytes);
    println!("checksum\t{}", stats.checksum);
}

fn one_path(args: Vec<String>, command: &str) -> Result<PathBuf> {
    if args.len() != 1 {
        return Err(microraptor::FastqError::Format(format!(
            "{command} requires exactly one path argument"
        )));
    }
    Ok(PathBuf::from(&args[0]))
}

fn required_path(path: Option<PathBuf>, command: &str) -> Result<PathBuf> {
    path.ok_or_else(|| microraptor::FastqError::Format(format!("{command} requires a path")))
}

fn required_value(iter: &mut impl Iterator<Item = String>, flag: &str) -> Result<String> {
    iter.next()
        .ok_or_else(|| microraptor::FastqError::Format(format!("{flag} requires a value")))
}

fn parse_u64_value(iter: &mut impl Iterator<Item = String>, flag: &str) -> Result<u64> {
    let value = required_value(iter, flag)?;
    value
        .parse::<u64>()
        .map_err(|_| microraptor::FastqError::Format(format!("{flag} requires an integer value")))
}

fn parse_usize_value(iter: &mut impl Iterator<Item = String>, flag: &str) -> Result<usize> {
    let value = required_value(iter, flag)?;
    value
        .parse::<usize>()
        .map_err(|_| microraptor::FastqError::Format(format!("{flag} requires an integer value")))
}

fn trim_line(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    line.strip_suffix(b"\r").unwrap_or(line)
}

fn print_help() {
    eprintln!(
        "microraptor <command>\n\ncommands:\n  stats [--format fastq|fasta] PATH\n  checksum --format fastq|fasta|sam [PATH|-]\n  fasta-index PATH\n  fasta-fetch PATH --fai PATH.fai --name REF --start N --end N\n  fasta-partitions PATH --fai PATH.fai --parts N --overlap N\n  fasta-chunks PATH --fai PATH.fai --name REF --start N --end N --chunk-bases N\n  verify-bgzf PATH"
    );
}
