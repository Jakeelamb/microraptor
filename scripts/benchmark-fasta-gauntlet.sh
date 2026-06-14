#!/usr/bin/env bash
set -euo pipefail

records="${MICRORAPTOR_FASTA_GAUNTLET_RECORDS:-100000}"
read_len="${MICRORAPTOR_FASTA_GAUNTLET_READ_LEN:-150}"
iters="${MICRORAPTOR_FASTA_GAUNTLET_ITERS:-3}"
workers="${MICRORAPTOR_BENCH_THREADS:-8}"
input_root="${MICRORAPTOR_FASTA_GAUNTLET_INPUT_DIR:-target/fasta-gauntlet-inputs}"
result_dir="${MICRORAPTOR_FASTA_GAUNTLET_RESULT_DIR:-target/bench-results/fasta-gauntlet}"
corpus_inputs="${MICRORAPTOR_FASTA_GAUNTLET_CORPUS_INPUTS:-}"
jsonl="${result_dir}/microraptor-fasta-gauntlet.jsonl"
external_tsv="${result_dir}/external-tools.tsv"
microraptor_memory_tsv="${result_dir}/microraptor-memory.tsv"
metadata="${result_dir}/metadata.md"
summary="${result_dir}/summary.md"

export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-${workers}}"
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-${workers}}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${workers}}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-${workers}}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-${workers}}"

mkdir -p "${input_root}" "${result_dir}"
: > "${jsonl}"
printf 'label\ttool\tstatus\telapsed_s\tmax_rss_kb\tcommand\n' > "${external_tsv}"
printf 'label\tstatus\telapsed_s\tmax_rss_kb\tcommand\n' > "${microraptor_memory_tsv}"

parse_corpus_inputs() {
  local value="$1"
  [[ -n "${value}" ]] || return 0
  if [[ "${value}" == *$'\n'* ]]; then
    printf '%s\n' "${value}"
  else
    # Backward-compatible space-separated form for simple paths.
    # shellcheck disable=SC2206
    local items=(${value})
    printf '%s\n' "${items[@]}"
  fi
}

command_version() {
  local command_name="$1"
  shift
  if command -v "${command_name}" >/dev/null 2>&1; then
    "$@" 2>&1 | sed -n '1,3p' || true
  else
    printf '%s not installed\n' "${command_name}"
  fi
}

cargo +nightly build --release --all-features --bin microraptor-bench --bin microraptor-fixture

make_fixture() {
  local label="$1"
  local fixture_records="$2"
  local fixture_read_len="$3"
  local pattern="$4"
  local layout="$5"
  local alphabet="$6"
  local out_dir="${input_root}/${label}"
  mkdir -p "${out_dir}"
  target/release/microraptor-fixture \
    --format fasta \
    --out-dir "${out_dir}" \
    --records "${fixture_records}" \
    --read-len "${fixture_read_len}" \
    --pattern "${pattern}" \
    --fasta-layout "${layout}" \
    --alphabet "${alphabet}"
}

make_fixture two-line-dna "${records}" "${read_len}" entropy two-line dna
make_fixture wrapped-dna "${records}" "${read_len}" entropy wrapped:60 dna
make_fixture many-tiny 100000 31 entropy two-line dna
make_fixture long-contigs 100 100000 entropy wrapped:80 dna
make_fixture protein 50000 300 entropy wrapped:80 protein

run_microraptor() {
  local label="$1"
  local path="$2"
  [[ -f "${path}" ]] || return 0
  printf 'running microraptor %s\n' "${label}"
  target/release/microraptor-bench \
    --format fasta \
    --mode parse \
    --input "${path}" \
    --iters "${iters}" \
    --workers "${workers}" \
    --json >> "${jsonl}"

  local status elapsed rss command_text
  command_text="$(printf 'target/release/microraptor-bench --format fasta --mode parse --input %q --iters 1 --workers %q --json' "${path}" "${workers}")"
  set +e
  /usr/bin/time -f '%e\t%M' -o "${result_dir}/.time.tmp" \
    target/release/microraptor-bench \
      --format fasta \
      --mode parse \
      --input "${path}" \
      --iters 1 \
      --workers "${workers}" \
      --json >/dev/null
  status="$?"
  set -e
  if [[ -s "${result_dir}/.time.tmp" ]]; then
    IFS=$'\t' read -r elapsed rss < "${result_dir}/.time.tmp"
  else
    elapsed=""
    rss=""
  fi
  if [[ "${status}" -eq 0 ]]; then
    printf '%s\tok\t%s\t%s\t%s\n' "${label}" "${elapsed}" "${rss}" "${command_text}" >> "${microraptor_memory_tsv}"
  else
    printf '%s\tfailed:%s\t%s\t%s\t%s\n' "${label}" "${status}" "${elapsed}" "${rss}" "${command_text}" >> "${microraptor_memory_tsv}"
  fi
}

run_external() {
  local label="$1"
  local tool="$2"
  shift 2
  if ! command -v "${tool}" >/dev/null 2>&1; then
    printf '%s\t%s\tskipped\t\t\t%s not installed\n' "${label}" "${tool}" "${tool}" >> "${external_tsv}"
    return 0
  fi

  local start end status elapsed rss command_text
  command_text="$(printf '%q ' "$@")"
  start="$(date +%s%N)"
  set +e
  /usr/bin/time -f '%e\t%M' -o "${result_dir}/.time.tmp" "$@" >/dev/null 2>"${result_dir}/.stderr.tmp"
  status="$?"
  set -e
  end="$(date +%s%N)"
  if [[ -s "${result_dir}/.time.tmp" ]]; then
    IFS=$'\t' read -r elapsed rss < "${result_dir}/.time.tmp"
  else
    elapsed="$(awk -v s="${start}" -v e="${end}" 'BEGIN { printf "%.6f", (e - s) / 1000000000.0 }')"
    rss=""
  fi
  if [[ "${status}" -eq 0 ]]; then
    printf '%s\t%s\tok\t%s\t%s\t%s\n' "${label}" "${tool}" "${elapsed}" "${rss}" "${command_text}" >> "${external_tsv}"
  else
    printf '%s\t%s\tfailed:%s\t%s\t%s\t%s\n' "${label}" "${tool}" "${status}" "${elapsed}" "${rss}" "${command_text}" >> "${external_tsv}"
  fi
}

run_external_samtools_faidx() {
  local label="$1"
  local path="$2"
  command -v samtools >/dev/null 2>&1 || {
    printf '%s\tsamtools\tskipped\t\t\tsamtools not installed\n' "${label}" >> "${external_tsv}"
    return 0
  }
  local tmp="${result_dir}/$(basename "${path}").faidx.tmp"
  cp "${path}" "${tmp}"
  run_external "${label}" samtools samtools faidx "${tmp}"
  rm -f "${tmp}" "${tmp}.fai" "${tmp}.gzi"
}

run_path_suite() {
  local label="$1"
  local path="$2"
  run_microraptor "${label}" "${path}"
  run_external "seqkit stats ${label}" seqkit seqkit stats "${path}"
  run_external "seqtk comp ${label}" seqtk seqtk comp "${path}"
  case "${path}" in
    *.gz | *.bgz)
      run_external "bgzip decompress ${label}" bgzip bgzip -dc "${path}"
      ;;
  esac
  case "${path}" in
    *.fasta | *.fa | *.fna | *.faa | *.bgz)
      run_external_samtools_faidx "samtools faidx ${label}" "${path}"
      ;;
  esac
}

for class_dir in "${input_root}"/*; do
  [[ -d "${class_dir}" ]] || continue
  class="$(basename "${class_dir}")"
  run_path_suite "${class}/raw" "${class_dir}/single.fasta"
  run_path_suite "${class}/gzip" "${class_dir}/single.fasta.gz"
  run_path_suite "${class}/bgzf" "${class_dir}/single.fasta.bgz"
done

while IFS= read -r corpus_input; do
  [[ -n "${corpus_input}" ]] || continue
  [[ -f "${corpus_input}" ]] || continue
  run_path_suite "corpus/$(basename "${corpus_input}")" "${corpus_input}"
done < <(parse_corpus_inputs "${corpus_inputs}")

{
  printf '# FASTA Gauntlet Metadata\n\n'
  printf -- '- generated_at_utc: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf -- '- git_commit: %s\n' "$(git rev-parse --short HEAD 2>/dev/null || printf unknown)"
  printf -- '- git_dirty: %s\n' "$(if [[ -n "$(git status --short 2>/dev/null)" ]]; then printf true; else printf false; fi)"
  printf -- '- rustc: %s\n' "$(rustc --version)"
  printf -- '- cargo: %s\n' "$(cargo --version)"
  printf -- '- records: %s\n' "${records}"
  printf -- '- read_len: %s\n' "${read_len}"
  printf -- '- iters: %s\n' "${iters}"
  printf -- '- thread_cap: %s\n' "${workers}"
  printf '\n## Tool Versions\n\n'
  printf '### seqkit\n\n```text\n'
  command_version seqkit seqkit version
  printf '```\n\n### seqtk\n\n```text\n'
  command_version seqtk seqtk
  printf '```\n\n### samtools\n\n```text\n'
  command_version samtools samtools --version
  printf '```\n\n### bgzip\n\n```text\n'
  command_version bgzip bgzip --version
  printf '```\n'
} > "${metadata}"

{
  printf '# FASTA Benchmark Gauntlet\n\n'
  printf 'Generated by `%s`.\n\n' "$0"
  printf 'This gauntlet covers synthetic two-line DNA, wrapped DNA, many tiny records, long wrapped contigs, protein FASTA, raw/gzip/BGZF transport, optional local corpus FASTA files, and command-line comparator timings when installed. Compression tools are third-party backends; microraptor rows measure parser/orchestration behavior over the selected transport.\n\n'
  printf 'Microraptor JSONL: [`microraptor-fasta-gauntlet.jsonl`](microraptor-fasta-gauntlet.jsonl)\n\n'
  printf 'Microraptor memory/RSS smoke rows: [`microraptor-memory.tsv`](microraptor-memory.tsv)\n\n'
  printf 'External tools: [`external-tools.tsv`](external-tools.tsv)\n\n'
  printf 'Metadata: [`metadata.md`](metadata.md)\n'
} > "${summary}"

rm -f "${result_dir}/.time.tmp" "${result_dir}/.stderr.tmp"

printf 'wrote %s\n' "${jsonl}"
printf 'wrote %s\n' "${external_tsv}"
printf 'wrote %s\n' "${microraptor_memory_tsv}"
printf 'wrote %s\n' "${metadata}"
printf 'wrote %s\n' "${summary}"
