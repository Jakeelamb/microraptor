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

run_microraptor_paired() {
  local label="$1"
  local first="$2"
  local second="$3"
  [[ -f "${first}" && -f "${second}" ]] || return 0

  printf 'running microraptor %s: %s %s\n' "${label}" "${first}" "${second}"
  target/release/microraptor-bench \
    --paired-inputs "${first}" "${second}" \
    --iters "${iters}" \
    --workers "${workers}" \
    --json >> "${jsonl}"

  {
    printf '### %s\n\n' "${label}"
    printf '```text\n'
    target/release/microraptor-bench \
      --paired-inputs "${first}" "${second}" \
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
      set +e
      /usr/bin/time -f 'elapsed_s\t%e' "$@" 2>&1
      local status="$?"
      set -e
      if [[ "${status}" -ne 0 ]]; then
        printf 'exit_status\t%s\n' "${status}"
      fi
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
run_microraptor_paired "paired/raw" "${input_dir}/r1.fastq" "${input_dir}/r2.fastq"
run_microraptor_paired "paired/gzip" "${input_dir}/r1.fastq.gz" "${input_dir}/r2.fastq.gz"
run_microraptor_paired "paired/bgzf" "${input_dir}/r1.fastq.bgz" "${input_dir}/r2.fastq.bgz"

{
  printf '## external tools\n\n'
} >> "${md}"

run_external "seqkit stats single/raw" seqkit seqkit stats "${input_dir}/single.fastq"
run_external "seqkit stats single/gzip" seqkit seqkit stats "${input_dir}/single.fastq.gz"
run_external "seqkit stats single/bgzf" seqkit seqkit stats "${input_dir}/single.fastq.bgz"
run_external "seqkit stats paired/r1/raw" seqkit seqkit stats "${input_dir}/r1.fastq"
run_external "seqkit stats paired/r2/raw" seqkit seqkit stats "${input_dir}/r2.fastq"
run_external "seqkit stats paired/r1/gzip" seqkit seqkit stats "${input_dir}/r1.fastq.gz"
run_external "seqkit stats paired/r2/gzip" seqkit seqkit stats "${input_dir}/r2.fastq.gz"
run_external "seqkit stats paired/r1/bgzf" seqkit seqkit stats "${input_dir}/r1.fastq.bgz"
run_external "seqkit stats paired/r2/bgzf" seqkit seqkit stats "${input_dir}/r2.fastq.bgz"
run_external "seqtk size single/raw" seqtk seqtk size "${input_dir}/single.fastq"
run_external "seqtk size single/gzip" seqtk seqtk size "${input_dir}/single.fastq.gz"
run_external "seqtk size single/bgzf" seqtk seqtk size "${input_dir}/single.fastq.bgz"
run_external "seqtk fqchk single/raw" seqtk seqtk fqchk "${input_dir}/single.fastq"
run_external "bgzip test single/bgzf" bgzip bgzip -t "${input_dir}/single.fastq.bgz"
run_external "bgzip decompress single/bgzf" bgzip bash -lc \
  "bgzip -dc '${input_dir}/single.fastq.bgz' >/dev/null"
run_external "samtools import single/raw" samtools samtools import \
  -0 "${input_dir}/single.fastq" \
  -o /dev/null \
  -O BAM \
  -@ "${workers}"
run_external "samtools import paired/raw" samtools samtools import \
  -1 "${input_dir}/r1.fastq" \
  -2 "${input_dir}/r2.fastq" \
  -o /dev/null \
  -O BAM \
  -@ "${workers}"
run_external "samtools import paired/gzip" samtools samtools import \
  -1 "${input_dir}/r1.fastq.gz" \
  -2 "${input_dir}/r2.fastq.gz" \
  -o /dev/null \
  -O BAM \
  -@ "${workers}"
run_external "samtools import paired/bgzf" samtools samtools import \
  -1 "${input_dir}/r1.fastq.bgz" \
  -2 "${input_dir}/r2.fastq.bgz" \
  -o /dev/null \
  -O BAM \
  -@ "${workers}"
run_external "fastp paired/raw" fastp bash -lc \
  "fastp --in1 '${input_dir}/r1.fastq' --in2 '${input_dir}/r2.fastq' --stdout --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering --thread '${workers}' --json '${result_dir}/fastp.json' --html '${result_dir}/fastp.html' >/dev/null"
run_external "fastp paired/gzip" fastp bash -lc \
  "fastp --in1 '${input_dir}/r1.fastq.gz' --in2 '${input_dir}/r2.fastq.gz' --stdout --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering --thread '${workers}' --json '${result_dir}/fastp-gzip.json' --html '${result_dir}/fastp-gzip.html' >/dev/null"
run_external "fastp paired/bgzf" fastp bash -lc \
  "fastp --in1 '${input_dir}/r1.fastq.bgz' --in2 '${input_dir}/r2.fastq.bgz' --stdout --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering --thread '${workers}' --json '${result_dir}/fastp-bgzf.json' --html '${result_dir}/fastp-bgzf.html' >/dev/null"

printf 'wrote %s\n' "${jsonl}"
printf 'wrote %s\n' "${md}"
