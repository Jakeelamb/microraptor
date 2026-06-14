#!/usr/bin/env bash
set -euo pipefail

project_dir="${MICRORAPTOR_FASTA_PEER_PROJECT_DIR:-target/fasta-peer-bench}"
out_dir="${MICRORAPTOR_FASTA_PEER_OUT_DIR:-docs/benchmarks/fasta-peers}"
iters="${MICRORAPTOR_FASTA_PEER_ITERS:-5}"
records="${MICRORAPTOR_FASTA_PEER_RECORDS:-100000}"
read_len="${MICRORAPTOR_FASTA_PEER_READ_LEN:-150}"
input="${MICRORAPTOR_FASTA_PEER_INPUT:-}"
raw_tsv="${project_dir}/fasta-library-peers.tsv"
tsv="${out_dir}/fasta-library-peers.tsv"
fig_dir="${out_dir}/figures"
peer_svg="${fig_dir}/fasta-library-peer-bases-throughput.svg"
summary="${out_dir}/summary.md"
metadata="${out_dir}/metadata.md"
microraptor_features="${MICRORAPTOR_FASTA_PEER_MICRORAPTOR_FEATURES:-}"
cargo_command="${MICRORAPTOR_FASTA_PEER_CARGO:-cargo}"
consumer="${MICRORAPTOR_FASTA_PEER_CONSUMER:-light}"
compression="${MICRORAPTOR_FASTA_PEER_COMPRESSION:-raw}"
bench_threads="${MICRORAPTOR_BENCH_THREADS:-8}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/benchmark-common.sh"
microraptor_set_thread_cap "${bench_threads}"

mkdir -p "${project_dir}/src" "${out_dir}" "${fig_dir}"

microraptor_dependency="$(microraptor_dependency_spec "${microraptor_features}")"

cat > "${project_dir}/Cargo.toml" <<EOF
[package]
name = "microraptor-fasta-peer-bench"
version = "0.0.0"
edition = "2024"
publish = false

[dependencies]
bio = "=4.0.0"
flate2 = "=1.1.9"
libdeflater = "=1.25.2"
${microraptor_dependency}
seq_io = "=0.3.4"
EOF

cat > "${project_dir}/src/main.rs" <<'EOF'
use std::env;
use std::fs;
use std::hint::black_box;
use std::io::{BufRead, BufReader, Cursor, Read, Write};
use std::mem::MaybeUninit;
use std::path::PathBuf;
use std::sync::OnceLock;
use std::time::{Duration, Instant};

use bio::io::fasta as bio_fasta;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression as GzipCompression;
use microraptor::benchutil::synthetic_fasta;
use microraptor::{
    count_two_line_fasta_bytes, count_two_line_fasta_read, visit_fasta_bytes,
    visit_two_line_fasta_bytes, visit_two_line_fasta_read, FastaConfig, FastaReader,
};

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

    fn observe(&mut self, seq: &[u8]) {
        self.records += 1;
        self.bases += seq.len() as u64;
        match consumer() {
            Consumer::FullChecksum => {
                self.checksum = mix_bytes(self.checksum, seq);
            }
            Consumer::LightAccounting => {
                self.checksum = mix_record_shape(self.checksum, seq);
            }
        }
    }
}

fn consumer() -> Consumer {
    static CONSUMER: OnceLock<Consumer> = OnceLock::new();
    *CONSUMER.get_or_init(|| match env::var("MICRORAPTOR_FASTA_PEER_CONSUMER") {
        Ok(value) if value == "light" => Consumer::LightAccounting,
        Ok(value) if value == "full" => Consumer::FullChecksum,
        Ok(value) if !value.is_empty() => {
            panic!("unsupported MICRORAPTOR_FASTA_PEER_CONSUMER={value}; expected full or light")
        }
        _ => Consumer::LightAccounting,
    })
}

fn fasta_batch_records() -> usize {
    static BATCH_RECORDS: OnceLock<usize> = OnceLock::new();
    *BATCH_RECORDS.get_or_init(|| match env::var("MICRORAPTOR_FASTA_PEER_BATCH_RECORDS") {
        Ok(value) if !value.is_empty() => value
            .parse()
            .unwrap_or_else(|_| panic!("invalid MICRORAPTOR_FASTA_PEER_BATCH_RECORDS={value}")),
        _ => 4096,
    })
}

fn microraptor_config() -> FastaConfig {
    FastaConfig {
        batch_records: fasta_batch_records(),
        ..FastaConfig::default()
    }
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
            let raw = synthetic_fasta(config.records, config.read_len);
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
            "microraptor-resident-visitor",
            &input,
            config.iters,
            config.compression,
            parse_microraptor_resident_visitor,
        )?,
        measure(
            "microraptor-two-line-resident-visitor",
            &input,
            config.iters,
            config.compression,
            parse_microraptor_two_line_resident_visitor,
        )?,
        measure(
            "microraptor-two-line-resident-counter",
            &input,
            config.iters,
            config.compression,
            parse_microraptor_two_line_resident_counter,
        )?,
        measure(
            "microraptor-two-line-stream-visitor",
            &input,
            config.iters,
            config.compression,
            parse_microraptor_two_line_stream_visitor,
        )?,
        measure(
            "microraptor-two-line-stream-counter",
            &input,
            config.iters,
            config.compression,
            parse_microraptor_two_line_stream_counter,
        )?,
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
    if config.compression == InputCompression::Gzip {
        rows.insert(
            1,
            measure(
                "microraptor-libdeflate-resident-visitor",
                &input,
                config.iters,
                config.compression,
                parse_microraptor_libdeflate_resident_visitor,
            )?,
        );
        rows.insert(
            2,
            measure(
                "microraptor-libdeflate-two-line-visitor",
                &input,
                config.iters,
                config.compression,
                parse_microraptor_libdeflate_two_line_visitor,
            )?,
        );
        rows.insert(
            3,
            measure(
                "microraptor-libdeflate-two-line-counter",
                &input,
                config.iters,
                config.compression,
                parse_microraptor_libdeflate_two_line_counter,
            )?,
        );
    }
    rows.extend([
        measure("seq_io", &input, config.iters, config.compression, parse_seq_io)?,
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
    let mut out = PathBuf::from("fasta-library-peers.tsv");
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
        "microraptor-fasta-peer-bench [--input PATH] [--records N] [--read-len N] [--iters N] [--out PATH] [--compression raw|gzip]"
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

fn input_read(input: &[u8], compression: InputCompression) -> Box<dyn Read + '_> {
    match compression {
        InputCompression::Raw => Box::new(Cursor::new(input)),
        InputCompression::Gzip => Box::new(MultiGzDecoder::new(Cursor::new(input))),
    }
}

fn gzip_bytes(input: &[u8]) -> AppResult<Vec<u8>> {
    let mut encoder = GzEncoder::new(Vec::new(), GzipCompression::default());
    encoder.write_all(input)?;
    Ok(encoder.finish()?)
}

fn parse_microraptor(input: &[u8], compression: InputCompression) -> AppResult<Stats> {
    let mut reader = FastaReader::with_config(input_reader(input, compression), microraptor_config());
    let mut stats = Stats::new();
    while let Some(batch) = reader.next_batch()? {
        for record in batch.records() {
            stats.observe(record.seq());
        }
    }
    Ok(stats)
}

fn parse_microraptor_visitor(input: &[u8], compression: InputCompression) -> AppResult<Stats> {
    let mut reader = FastaReader::with_config(input_reader(input, compression), microraptor_config());
    let mut stats = Stats::new();
    reader.visit_records(|record| {
        stats.observe(record.seq());
        Ok(())
    })?;
    Ok(stats)
}

fn parse_microraptor_resident_visitor(input: &[u8], compression: InputCompression) -> AppResult<Stats> {
    let decoded;
    let bytes = match compression {
        InputCompression::Raw => input,
        InputCompression::Gzip => {
            let mut decoder = MultiGzDecoder::new(Cursor::new(input));
            decoded = {
                let mut decoded = Vec::with_capacity(initial_gzip_output_capacity(input));
                decoder.read_to_end(&mut decoded)?;
                decoded
            };
            &decoded
        }
    };
    let mut stats = Stats::new();
    visit_fasta_bytes(bytes, |record| {
        stats.observe(record.seq());
        Ok(())
    })?;
    Ok(stats)
}

fn parse_microraptor_two_line_resident_visitor(
    input: &[u8],
    compression: InputCompression,
) -> AppResult<Stats> {
    let decoded;
    let bytes = match compression {
        InputCompression::Raw => input,
        InputCompression::Gzip => {
            let mut decoder = MultiGzDecoder::new(Cursor::new(input));
            decoded = {
                let mut decoded = Vec::with_capacity(initial_gzip_output_capacity(input));
                decoder.read_to_end(&mut decoded)?;
                decoded
            };
            &decoded
        }
    };
    let mut stats = Stats::new();
    visit_two_line_fasta_bytes(bytes, |record| {
        stats.observe(record.seq());
        Ok(())
    })?;
    Ok(stats)
}

fn parse_microraptor_two_line_resident_counter(
    input: &[u8],
    compression: InputCompression,
) -> AppResult<Stats> {
    let decoded;
    let bytes = match compression {
        InputCompression::Raw => input,
        InputCompression::Gzip => {
            let mut decoder = MultiGzDecoder::new(Cursor::new(input));
            decoded = {
                let mut decoded = Vec::with_capacity(initial_gzip_output_capacity(input));
                decoder.read_to_end(&mut decoded)?;
                decoded
            };
            &decoded
        }
    };
    let stats = count_two_line_fasta_bytes(bytes)?;
    Ok(Stats {
        records: stats.records,
        bases: stats.bases,
        checksum: stats.checksum,
    })
}

fn parse_microraptor_two_line_stream_visitor(
    input: &[u8],
    compression: InputCompression,
) -> AppResult<Stats> {
    let mut stats = Stats::new();
    visit_two_line_fasta_read(input_read(input, compression), |record| {
        stats.observe(record.seq());
        Ok(())
    })?;
    Ok(stats)
}

fn parse_microraptor_two_line_stream_counter(
    input: &[u8],
    compression: InputCompression,
) -> AppResult<Stats> {
    let stats = count_two_line_fasta_read(input_read(input, compression))?;
    Ok(Stats {
        records: stats.records,
        bases: stats.bases,
        checksum: stats.checksum,
    })
}

fn parse_microraptor_libdeflate_resident_visitor(
    input: &[u8],
    compression: InputCompression,
) -> AppResult<Stats> {
    let decoded = match compression {
        InputCompression::Raw => input.to_vec(),
        InputCompression::Gzip => decompress_gzip_libdeflate(input)?,
    };
    let mut stats = Stats::new();
    visit_fasta_bytes(&decoded, |record| {
        stats.observe(record.seq());
        Ok(())
    })?;
    Ok(stats)
}

fn parse_microraptor_libdeflate_two_line_visitor(
    input: &[u8],
    compression: InputCompression,
) -> AppResult<Stats> {
    let decoded = match compression {
        InputCompression::Raw => input.to_vec(),
        InputCompression::Gzip => decompress_gzip_libdeflate(input)?,
    };
    let mut stats = Stats::new();
    visit_two_line_fasta_bytes(&decoded, |record| {
        stats.observe(record.seq());
        Ok(())
    })?;
    Ok(stats)
}

fn parse_microraptor_libdeflate_two_line_counter(
    input: &[u8],
    compression: InputCompression,
) -> AppResult<Stats> {
    let decoded = match compression {
        InputCompression::Raw => input.to_vec(),
        InputCompression::Gzip => decompress_gzip_libdeflate(input)?,
    };
    let stats = count_two_line_fasta_bytes(&decoded)?;
    Ok(Stats {
        records: stats.records,
        bases: stats.bases,
        checksum: stats.checksum,
    })
}

fn parse_seq_io(input: &[u8], compression: InputCompression) -> AppResult<Stats> {
    let mut reader = seq_io::fasta::Reader::new(input_reader(input, compression));
    let mut stats = Stats::new();
    while let Some(record) = reader.next() {
        let record = record?;
        let seq = record.full_seq();
        stats.observe(seq.as_ref());
    }
    Ok(stats)
}

fn parse_bio(input: &[u8], compression: InputCompression) -> AppResult<Stats> {
    let reader = bio_fasta::Reader::new(input_reader(input, compression));
    let mut stats = Stats::new();
    for record in reader.records() {
        let record = record?;
        stats.observe(record.seq());
    }
    Ok(stats)
}

fn decompress_gzip_libdeflate(input: &[u8]) -> AppResult<Vec<u8>> {
    let mut out = Vec::<MaybeUninit<u8>>::with_capacity(initial_gzip_output_capacity(input));
    let mut decompressor = libdeflater::Decompressor::new();
    loop {
        let capacity = out.capacity().max(1);
        let output = unsafe {
            std::slice::from_raw_parts_mut(out.as_mut_ptr().cast::<u8>(), capacity)
        };
        match decompressor.gzip_decompress(input, output) {
            Ok(n) => {
                let ptr = out.as_mut_ptr().cast::<u8>();
                let cap = out.capacity();
                std::mem::forget(out);
                return Ok(unsafe { Vec::from_raw_parts(ptr, n, cap) });
            }
            Err(libdeflater::DecompressionError::InsufficientSpace) => {
                let next = out
                    .capacity()
                    .checked_mul(2)
                    .ok_or("gzip output size exceeds usize range")?;
                out.reserve_exact(next.saturating_sub(out.capacity()).max(1));
            }
            Err(err) => return Err(format!("libdeflate gzip inflate failed: {err}").into()),
        }
    }
}

fn initial_gzip_output_capacity(input: &[u8]) -> usize {
    let isize = input
        .len()
        .checked_sub(4)
        .and_then(|start| input.get(start..))
        .map(|tail| u32::from_le_bytes([tail[0], tail[1], tail[2], tail[3]]) as usize)
        .unwrap_or(0);
    isize.max(input.len().saturating_mul(2)).max(1024)
}

fn mix_bytes(mut state: u64, bytes: &[u8]) -> u64 {
    for &byte in bytes {
        state ^= byte as u64;
        state = state.wrapping_mul(0x0000_0100_0000_01b3);
    }
    state
}

fn mix_record_shape(mut state: u64, seq: &[u8]) -> u64 {
    state ^= seq.len() as u64;
    state = state.rotate_left(5).wrapping_mul(0x0000_0100_0000_01b3);
    state ^= seq.first().copied().unwrap_or_default() as u64;
    state ^= (seq.last().copied().unwrap_or_default() as u64) << 8;
    state
}
EOF

args=(--iters "${iters}" --out "${raw_tsv}" --compression "${compression}")
if [[ -n "${input}" ]]; then
  args+=(--input "${input}")
else
  args+=(--records "${records}" --read-len "${read_len}")
fi

read -r -a cargo_cmd <<< "${cargo_command}"
MICRORAPTOR_FASTA_PEER_CONSUMER="${consumer}" "${cargo_cmd[@]}" run --release --manifest-path "${project_dir}/Cargo.toml" -- "${args[@]}"

microraptor_write_sanitized_file "${raw_tsv}" "${tsv}"

awk -F '\t' '
  NR > 1 {
    tool[++n] = $1
    bases_s[n] = $6 + 0
    if (bases_s[n] > max) {
      max = bases_s[n]
    }
  }
  END {
    width = 880
    left = 190
    bar_max = 500
    row_h = 34
    top = 58
    height = top + n * row_h + 44
    print "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" width "\" height=\"" height "\" viewBox=\"0 0 " width " " height "\">"
    print "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>"
    print "<text x=\"24\" y=\"28\" font-family=\"Arial, sans-serif\" font-size=\"18\" font-weight=\"700\">Rust FASTA parser peer throughput</text>"
    print "<text x=\"24\" y=\"46\" font-family=\"Arial, sans-serif\" font-size=\"12\" fill=\"#555\">bases/s over the same in-memory FASTA byte buffer</text>"
    for (i = 1; i <= n; i++) {
      y = top + (i - 1) * row_h
      bar = max > 0 ? int((bases_s[i] / max) * bar_max) : 0
      printf "<text x=\"24\" y=\"%d\" font-family=\"Arial, sans-serif\" font-size=\"12\" fill=\"#222\">%s</text>\n", y + 18, tool[i]
      printf "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"20\" fill=\"#23645b\"/>\n", left, y + 3, bar
      printf "<text x=\"%d\" y=\"%d\" font-family=\"Arial, sans-serif\" font-size=\"12\" fill=\"#222\">%.2f Gbases/s</text>\n", left + bar + 8, y + 18, bases_s[i] / 1000000000.0
    }
    print "</svg>"
  }
' "${tsv}" > "${peer_svg}"

{
  printf '# FASTA Rust Library Peer Benchmark\n\n'
  printf 'Generated by `%s`.\n\n' "$0"
  if [[ -n "${input}" ]]; then
    printf 'Input: `%s`\n\n' "$(microraptor_display_path "${input}")"
  else
    printf 'Input: deterministic synthetic FASTA, `%s` records, sequence length `%s`.\n\n' "${records}" "${read_len}"
  fi
  if [[ "${compression}" == "gzip" ]]; then
    printf 'The benchmark reads one gzip-compressed in-memory FASTA byte buffer through each Rust parser. Streaming rows use the same flate2 decoder directly. Resident rows first decompress to a byte buffer and then use microraptor resident visitors. Rows containing `two-line` use the strict `>header`/`sequence` fast path for canonical two-line FASTA; use the non-two-line resident row for ordinary multiline FASTA. The `counter` row is the light-accounting scan path used for count/bases/checksum workloads.\n\n'
  else
    printf 'The benchmark reads one in-memory raw FASTA byte buffer through each Rust parser. It is parser-library evidence only and does not compare indexing, filtering, or command-line workflow behavior. Rows containing `two-line` use the strict `>header`/`sequence` fast path for canonical two-line FASTA; use the non-two-line resident row for ordinary multiline FASTA. The `counter` row is the light-accounting scan path used for count/bases/checksum workloads.\n\n'
  fi
  printf '| tool | records | bases | best ms | records/s | bases/s | checksum |\n'
  printf '| --- | ---: | ---: | ---: | ---: | ---: | ---: |\n'
  awk -F '\t' 'NR > 1 {
    printf "| `%s` | %s | %s | %.3f | %s | %s | `%s` |\n", $1, $2, $3, $4, $5, $6, $7
  }' "${tsv}"
  printf '\nFigure: [`figures/fasta-library-peer-bases-throughput.svg`](figures/fasta-library-peer-bases-throughput.svg)\n\n'
  printf 'Metadata: [`metadata.md`](metadata.md)\n'
} > "${summary}"

{
  printf '# FASTA Rust Library Peer Benchmark Metadata\n\n'
  printf -- '- generated_at_utc: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf -- '- rustc: %s\n' "$(rustc --version)"
  printf -- '- cargo: %s\n' "$(cargo --version)"
  printf -- '- iterations: %s\n' "${iters}"
  printf -- '- consumer: %s\n' "${consumer}"
  printf -- '- compression: %s\n' "${compression}"
  printf -- '- thread_cap: %s\n' "${bench_threads}"
  printf -- '- microraptor_features: %s\n' "${microraptor_features:-default}"
  printf -- '- microraptor_batch_records: %s\n' "${MICRORAPTOR_FASTA_PEER_BATCH_RECORDS:-4096}"
  if [[ -n "${input}" ]]; then
    printf -- '- input: %s\n' "$(microraptor_display_path "${input}")"
  else
    printf -- '- input: synthetic:%sx%s\n' "${records}" "${read_len}"
  fi
  printf '\n## Resolved Crate Versions\n\n'
  awk '
    $1 == "name" {
      name = $3
      gsub(/"/, "", name)
      keep = name == "microraptor" || name == "seq_io" || name == "bio" || name == "libdeflater"
    }
    keep && $1 == "version" {
      version = $3
      gsub(/"/, "", version)
      printf "- %s: %s\n", name, version
      keep = 0
    }
  ' "${project_dir}/Cargo.lock"
} > "${metadata}"

printf 'wrote %s\n' "${tsv}"
printf 'wrote %s\n' "${summary}"
printf 'wrote %s\n' "${peer_svg}"
printf 'wrote %s\n' "${metadata}"
