#!/usr/bin/env bash
set -euo pipefail

records="${MICRORAPTOR_GAUNTLET_RECORDS:-100000}"
read_len="${MICRORAPTOR_GAUNTLET_READ_LEN:-150}"
iters="${MICRORAPTOR_GAUNTLET_ITERS:-3}"
workers="${MICRORAPTOR_WORKERS:-$(nproc)}"
input_dir="${MICRORAPTOR_GAUNTLET_INPUT_DIR:-target/bench-inputs}"
result_dir="${MICRORAPTOR_GAUNTLET_RESULT_DIR:-target/bench-results}"

mkdir -p "${input_dir}" "${result_dir}"

cargo build --release --all-features --bin microraptor-bench --bin microraptor-fixture
target/release/microraptor-fixture \
  --out-dir "${input_dir}" \
  --records "${records}" \
  --read-len "${read_len}"

jsonl="${result_dir}/microraptor-gauntlet.jsonl"
md="${result_dir}/microraptor-gauntlet.md"
: > "${jsonl}"

{
  printf '# microraptor benchmark gauntlet\n\n'
  printf -- '- records: %s\n' "${records}"
  printf -- '- read_len: %s\n' "${read_len}"
  printf -- '- iters: %s\n' "${iters}"
  printf -- '- workers: %s\n\n' "${workers}"
  printf '## microraptor\n\n'
} > "${md}"

run_microraptor() {
  local label="$1"
  local path="$2"
  [[ -f "${path}" ]] || return 0

  printf 'running microraptor %s: %s\n' "${label}" "${path}"
  target/release/microraptor-bench \
    --input "${path}" \
    --iters "${iters}" \
    --workers "${workers}" \
    --json >> "${jsonl}"

  {
    printf '### %s\n\n' "${label}"
    printf '```text\n'
    target/release/microraptor-bench \
      --input "${path}" \
      --iters 1 \
      --workers "${workers}"
    printf '```\n\n'
  } >> "${md}"
}

run_external() {
  local label="$1"
  local command_name="$2"
  shift 2

  {
    printf '### %s\n\n' "${label}"
    if command -v "${command_name}" >/dev/null 2>&1; then
      printf '```text\n'
      /usr/bin/time -f 'elapsed_s\t%e' "$@" 2>&1
      printf '```\n\n'
    else
      printf '`%s` not installed; skipped.\n\n' "${command_name}"
    fi
  } >> "${md}"
}

run_microraptor "single/raw" "${input_dir}/single.fastq"
run_microraptor "single/gzip" "${input_dir}/single.fastq.gz"
run_microraptor "single/bgzf" "${input_dir}/single.fastq.bgz"
run_microraptor "interleaved/raw" "${input_dir}/interleaved.fastq"
run_microraptor "interleaved/gzip" "${input_dir}/interleaved.fastq.gz"
run_microraptor "interleaved/bgzf" "${input_dir}/interleaved.fastq.bgz"
run_microraptor "paired/r1/raw" "${input_dir}/r1.fastq"
run_microraptor "paired/r2/raw" "${input_dir}/r2.fastq"

{
  printf '## external tools\n\n'
} >> "${md}"

run_external "seqkit stats single/raw" seqkit seqkit stats "${input_dir}/single.fastq"
run_external "seqkit stats single/gzip" seqkit seqkit stats "${input_dir}/single.fastq.gz"
run_external "fastp paired/raw" fastp fastp \
  --in1 "${input_dir}/r1.fastq" \
  --in2 "${input_dir}/r2.fastq" \
  --stdout \
  --disable_adapter_trimming \
  --disable_quality_filtering \
  --disable_length_filtering \
  --thread "${workers}"

printf 'wrote %s\n' "${jsonl}"
printf 'wrote %s\n' "${md}"
