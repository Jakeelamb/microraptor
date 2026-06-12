#!/usr/bin/env bash
set -euo pipefail

project_dir="${MICRORAPTOR_RUST_PEER_PROJECT_DIR:-target/rust-peer-bench}"
out_dir="${MICRORAPTOR_RUST_PEER_OUT_DIR:-docs/benchmarks/rust-peers}"
iters="${MICRORAPTOR_RUST_PEER_ITERS:-5}"
records="${MICRORAPTOR_RUST_PEER_RECORDS:-100000}"
read_len="${MICRORAPTOR_RUST_PEER_READ_LEN:-150}"
input="${MICRORAPTOR_RUST_PEER_INPUT:-}"
raw_tsv="${project_dir}/rust-library-peers.tsv"
tsv="${out_dir}/rust-library-peers.tsv"
fig_dir="${out_dir}/figures"
peer_svg="${fig_dir}/rust-library-peer-bases-throughput.svg"
summary="${out_dir}/summary.md"
metadata="${out_dir}/metadata.md"
microraptor_features="${MICRORAPTOR_RUST_PEER_MICRORAPTOR_FEATURES:-}"
cargo_command="${MICRORAPTOR_RUST_PEER_CARGO:-cargo}"

display_path() {
  if [[ -n "${HOME:-}" ]]; then
    awk -v home="${HOME}" 'BEGIN { value = ARGV[1]; ARGV[1] = ""; gsub(home, "~", value); print value }' "$1"
  else
    printf '%s\n' "$1"
  fi
}

mkdir -p "${project_dir}/src" "${out_dir}" "${fig_dir}"

if [[ -n "${microraptor_features}" ]]; then
  feature_list="$(printf '%s\n' "${microraptor_features}" | awk -F ',' '{
    for (i = 1; i <= NF; i++) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $i)
      if ($i != "") {
        printf "%s\"%s\"", sep, $i
        sep = ", "
      }
    }
  }')"
  microraptor_dependency="microraptor = { path = \"../..\", features = [${feature_list}] }"
else
  microraptor_dependency='microraptor = { path = "../.." }'
fi

cat > "${project_dir}/Cargo.toml" <<EOF
[package]
name = "microraptor-rust-peer-bench"
version = "0.0.0"
edition = "2024"
publish = false

[dependencies]
bio = "=4.0.0"
${microraptor_dependency}
noodles-fastq = "=0.23.0"
seq_io = "=0.3.4"
EOF

cat > "${project_dir}/src/main.rs" <<'EOF'
use std::env;
use std::fs;
use std::hint::black_box;
use std::io::Cursor;
use std::path::PathBuf;
use std::time::{Duration, Instant};

use bio::io::fastq as bio_fastq;
use microraptor::benchutil::synthetic_fastq;
use microraptor::{FastqConfig, FastqReader};
use noodles_fastq as noodles_fastq;
use seq_io::fastq::Record as SeqIoRecord;

type AppResult<T> = Result<T, Box<dyn std::error::Error>>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Stats {
    records: u64,
    bases: u64,
    checksum: u64,
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
        self.checksum = mix_bytes(self.checksum, seq);
        self.checksum = mix_bytes(self.checksum, qual);
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
}

fn main() -> AppResult<()> {
    let config = parse_args()?;
    let (source, input) = match config.input.as_ref() {
        Some(path) => (path.display().to_string(), fs::read(path)?),
        None => (
            format!("synthetic:{}x{}", config.records, config.read_len),
            synthetic_fastq(config.records, config.read_len),
        ),
    };

    let mut rows = vec![measure(
        "microraptor",
        &input,
        config.iters,
        parse_microraptor,
    )?];
    if env::var_os("MICRORAPTOR_RUST_PEER_DIAGNOSTICS").is_some() {
        rows.push(measure(
            "microraptor-no-validate",
            &input,
            config.iters,
            parse_microraptor_no_validate,
        )?);
        rows.push(measure(
            "microraptor-record-refs",
            &input,
            config.iters,
            parse_microraptor_record_refs,
        )?);
    }
    rows.extend([
        measure("seq_io", &input, config.iters, parse_seq_io)?,
        measure("noodles-fastq", &input, config.iters, parse_noodles_fastq)?,
        measure("bio", &input, config.iters, parse_bio)?,
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
        "tool\trecords\tbases\tbest_ms\trecords_s\tbases_s\tchecksum\tinput_bytes\titers\tsource\n",
    );
    for row in rows {
        let ns = row.best.as_nanos().max(1);
        let records_s = row.stats.records as u128 * 1_000_000_000 / ns;
        let bases_s = row.stats.bases as u128 * 1_000_000_000 / ns;
        out.push_str(&format!(
            "{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            row.tool,
            row.stats.records,
            row.stats.bases,
            row.best.as_secs_f64() * 1000.0,
            records_s,
            bases_s,
            row.stats.checksum,
            input.len(),
            config.iters,
            source
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

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--input" => input = Some(PathBuf::from(required_value(&mut args, "--input")?)),
            "--records" => records = required_value(&mut args, "--records")?.parse()?,
            "--read-len" => read_len = required_value(&mut args, "--read-len")?.parse()?,
            "--iters" => iters = required_value(&mut args, "--iters")?.parse()?,
            "--out" => out = PathBuf::from(required_value(&mut args, "--out")?),
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
    })
}

fn required_value(args: &mut impl Iterator<Item = String>, flag: &str) -> AppResult<String> {
    args.next()
        .ok_or_else(|| format!("missing value for {flag}").into())
}

fn print_help() {
    println!(
        "microraptor-rust-peer-bench [--input PATH] [--records N] [--read-len N] [--iters N] [--out PATH]"
    );
}

fn measure(
    tool: &'static str,
    input: &[u8],
    iters: usize,
    f: fn(&[u8]) -> AppResult<Stats>,
) -> AppResult<Row> {
    let mut best = Duration::MAX;
    let mut stats = None;
    for _ in 0..iters {
        let start = Instant::now();
        let run_stats = f(input)?;
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

fn parse_microraptor(input: &[u8]) -> AppResult<Stats> {
    let mut reader = FastqReader::with_config(
        Cursor::new(input),
        FastqConfig {
            validate: true,
            ..FastqConfig::default()
        },
    );
    let mut stats = Stats::new();
    while let Some(batch) = reader.next_batch()? {
        for record in batch.records() {
            stats.observe(record.seq(), record.qual());
        }
    }
    Ok(stats)
}

fn parse_microraptor_no_validate(input: &[u8]) -> AppResult<Stats> {
    let mut reader = FastqReader::with_config(
        Cursor::new(input),
        FastqConfig {
            validate: false,
            ..FastqConfig::default()
        },
    );
    let mut stats = Stats::new();
    while let Some(batch) = reader.next_batch()? {
        for record in batch.records() {
            stats.observe(record.seq(), record.qual());
        }
    }
    Ok(stats)
}

fn parse_microraptor_record_refs(input: &[u8]) -> AppResult<Stats> {
    let mut reader = FastqReader::with_config(
        Cursor::new(input),
        FastqConfig {
            validate: false,
            ..FastqConfig::default()
        },
    );
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

fn parse_seq_io(input: &[u8]) -> AppResult<Stats> {
    let mut reader = seq_io::fastq::Reader::new(Cursor::new(input));
    let mut stats = Stats::new();
    while let Some(record) = reader.next() {
        let record = record?;
        stats.observe(record.seq(), record.qual());
    }
    Ok(stats)
}

fn parse_noodles_fastq(input: &[u8]) -> AppResult<Stats> {
    let mut reader = noodles_fastq::io::Reader::new(Cursor::new(input));
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

fn parse_bio(input: &[u8]) -> AppResult<Stats> {
    let reader = bio_fastq::Reader::new(Cursor::new(input));
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
EOF

args=(--iters "${iters}" --out "${raw_tsv}")
if [[ -n "${input}" ]]; then
  args+=(--input "${input}")
else
  args+=(--records "${records}" --read-len "${read_len}")
fi

read -r -a cargo_cmd <<< "${cargo_command}"
"${cargo_cmd[@]}" run --release --manifest-path "${project_dir}/Cargo.toml" -- "${args[@]}"

if [[ -n "${HOME:-}" ]]; then
  awk -v home="${HOME}" '{ gsub(home, "~"); print }' "${raw_tsv}" > "${tsv}"
else
  cp "${raw_tsv}" "${tsv}"
fi

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
    left = 170
    bar_max = 520
    row_h = 34
    top = 58
    height = top + n * row_h + 44
    print "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" width "\" height=\"" height "\" viewBox=\"0 0 " width " " height "\">"
    print "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>"
    print "<text x=\"24\" y=\"28\" font-family=\"Arial, sans-serif\" font-size=\"18\" font-weight=\"700\">Rust FASTQ parser peer throughput</text>"
    print "<text x=\"24\" y=\"46\" font-family=\"Arial, sans-serif\" font-size=\"12\" fill=\"#555\">bases/s over the same in-memory raw FASTQ byte buffer</text>"
    for (i = 1; i <= n; i++) {
      y = top + (i - 1) * row_h
      bar = max > 0 ? int((bases_s[i] / max) * bar_max) : 0
      printf "<text x=\"24\" y=\"%d\" font-family=\"Arial, sans-serif\" font-size=\"12\" fill=\"#222\">%s</text>\n", y + 18, tool[i]
      printf "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"20\" fill=\"#315f8c\"/>\n", left, y + 3, bar
      printf "<text x=\"%d\" y=\"%d\" font-family=\"Arial, sans-serif\" font-size=\"12\" fill=\"#222\">%.2f Gbases/s</text>\n", left + bar + 8, y + 18, bases_s[i] / 1000000000.0
    }
    print "</svg>"
  }
' "${tsv}" > "${peer_svg}"

{
  printf '# Rust Library Peer Benchmark\n\n'
  printf 'Generated by `%s`.\n\n' "$0"
  if [[ -n "${input}" ]]; then
    printf 'Input: `%s`\n\n' "$(display_path "${input}")"
  else
    printf 'Input: deterministic synthetic FASTQ, `%s` records, read length `%s`.\n\n' "${records}" "${read_len}"
  fi
  printf 'The benchmark reads one in-memory raw FASTQ byte buffer through each Rust parser. It is parser-library evidence only: it does not compare gzip, BGZF, trimming, filtering, or command-line workflow behavior.\n\n'
  printf '| tool | records | bases | best ms | records/s | bases/s | checksum |\n'
  printf '| --- | ---: | ---: | ---: | ---: | ---: | ---: |\n'
  awk -F '\t' 'NR > 1 {
    printf "| `%s` | %s | %s | %.3f | %s | %s | `%s` |\n", $1, $2, $3, $4, $5, $6, $7
  }' "${tsv}"
  printf '\nFigure: [`figures/rust-library-peer-bases-throughput.svg`](figures/rust-library-peer-bases-throughput.svg)\n\n'
  printf 'Metadata: [`metadata.md`](metadata.md)\n'
} > "${summary}"

{
  printf '# Rust Library Peer Benchmark Metadata\n\n'
  printf -- '- generated_at_utc: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf -- '- rustc: %s\n' "$(rustc --version)"
  printf -- '- cargo: %s\n' "$(cargo --version)"
  printf -- '- iterations: %s\n' "${iters}"
  printf -- '- microraptor_features: %s\n' "${microraptor_features:-default}"
  if [[ -n "${input}" ]]; then
    printf -- '- input: %s\n' "$(display_path "${input}")"
  else
    printf -- '- input: synthetic:%sx%s\n' "${records}" "${read_len}"
  fi
  printf '\n## Resolved Crate Versions\n\n'
  awk '
    $1 == "name" {
      name = $3
      gsub(/"/, "", name)
      keep = name == "microraptor" || name == "seq_io" || name == "noodles-fastq" || name == "bio"
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
