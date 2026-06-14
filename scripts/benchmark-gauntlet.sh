#!/usr/bin/env bash
set -euo pipefail

records="${MICRORAPTOR_GAUNTLET_RECORDS:-100000}"
read_len="${MICRORAPTOR_GAUNTLET_READ_LEN:-150}"
iters="${MICRORAPTOR_GAUNTLET_ITERS:-3}"
workers="${MICRORAPTOR_WORKERS:-$(nproc)}"
input_dir="${MICRORAPTOR_GAUNTLET_INPUT_DIR:-target/bench-inputs}"
result_dir="${MICRORAPTOR_GAUNTLET_RESULT_DIR:-target/bench-results}"
corpus_inputs="${MICRORAPTOR_GAUNTLET_CORPUS_INPUTS:-}"
corpus_paired_inputs="${MICRORAPTOR_GAUNTLET_CORPUS_PAIRED_INPUTS:-}"
corpus_input_list=()
corpus_paired_input_list=()

parse_corpus_list() {
  local value="$1"
  local list_name="$2"
  local -n list_ref="${list_name}"
  list_ref=()

  [[ -n "${value}" ]] || return 0
  if [[ "${value}" == *$'\n'* ]]; then
    local item
    while IFS= read -r item; do
      [[ -n "${item}" ]] || continue
      list_ref+=("${item}")
    done <<< "${value}"
  else
    # Backward-compatible parser for the older space-separated env format.
    # New generated env files are newline-delimited so paths may contain spaces.
    # shellcheck disable=SC2206
    list_ref=(${value})
  fi
}

parse_corpus_list "${corpus_inputs}" corpus_input_list
parse_corpus_list "${corpus_paired_inputs}" corpus_paired_input_list

mkdir -p "${input_dir}" "${result_dir}"

cargo build --release --all-features --bin microraptor --bin microraptor-bench --bin microraptor-fixture
target/release/microraptor-fixture \
  --out-dir "${input_dir}" \
  --records "${records}" \
  --read-len "${read_len}"

jsonl="${result_dir}/microraptor-gauntlet.jsonl"
md="${result_dir}/microraptor-gauntlet.md"
metadata="${result_dir}/microraptor-gauntlet-metadata.md"
external_tsv="${result_dir}/external-tools.tsv"
external_parity_tsv="${result_dir}/external-parity.tsv"
: > "${jsonl}"
: > "${external_tsv}"
: > "${external_parity_tsv}"

command_version() {
  local command_name="$1"
  shift
  if command -v "${command_name}" >/dev/null 2>&1; then
    "$@" 2>&1 | sed -n '1,3p' || true
  else
    printf '%s not installed\n' "${command_name}"
  fi
}

{
  printf '# microraptor benchmark metadata\n\n'
  printf -- '- generated_at_utc: %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  printf -- '- git_commit: %s\n' "$(git rev-parse --short HEAD 2>/dev/null || printf unknown)"
  printf -- '- git_dirty: %s\n' "$(if [[ -n "$(git status --short 2>/dev/null)" ]]; then printf true; else printf false; fi)"
  printf -- '- rustc: %s\n' "$(rustc --version)"
  printf -- '- cargo: %s\n' "$(cargo --version)"
  printf -- '- cargo_nightly: %s\n' "$(cargo +nightly --version)"
  printf -- '- build_profile: release\n'
  printf -- '- feature_flags: all-features\n'
  printf -- '- kernel: %s\n' "$(uname -srmo)"
  printf -- '- cpu: %s\n' "$(lscpu | sed -n 's/^Model name:[[:space:]]*//p' | head -n 1)"
  printf -- '- logical_cpus: %s\n' "$(nproc)"
  printf -- '- memory: %s\n' "$(free -h | awk '/^Mem:/ { print $2 }')"
  printf -- '- filesystem: %s\n' "$(df -T . | awk 'NR == 2 { printf "%s", $2 }')"
  printf -- '- storage_available: %s\n' "$(df -h . | awk 'NR == 2 { printf "%s", $4 }')"
  printf -- '- records: %s\n' "${records}"
  printf -- '- read_len: %s\n' "${read_len}"
  printf -- '- iters: %s\n' "${iters}"
  printf -- '- workers: %s\n\n' "${workers}"
  printf '## Tool versions\n\n'
  printf '### microraptor\n\n```text\n'
  target/release/microraptor-bench --help 2>&1 | sed -n '1p' || true
  printf '```\n\n'
  printf '### seqkit\n\n```text\n'
  command_version seqkit seqkit version
  printf '```\n\n'
  printf '### seqtk\n\n```text\n'
  command_version seqtk seqtk
  printf '```\n\n'
  printf '### samtools\n\n```text\n'
  command_version samtools samtools --version
  printf '```\n\n'
  printf '### bgzip\n\n```text\n'
  command_version bgzip bgzip --version
  printf '```\n\n'
  printf '### fastp\n\n```text\n'
  command_version fastp fastp --version
  printf '```\n'
} > "${metadata}"

{
  printf '# microraptor benchmark gauntlet\n\n'
  printf -- '- records: %s\n' "${records}"
  printf -- '- read_len: %s\n' "${read_len}"
  printf -- '- iters: %s\n' "${iters}"
  printf -- '- workers: %s\n\n' "${workers}"
  printf '## microraptor\n\n'
} > "${md}"

printf 'label\ttool\tstatus\telapsed_s\tcommand\n' > "${external_tsv}"
printf 'label\ttool\tparity_status\texpected_records\texpected_bases\texpected_checksum\tobserved_records\tobserved_bases\tobserved_checksum\tnotes\n' > "${external_parity_tsv}"

record_external_parity_timing_only() {
  local label="$1"
  local command_name="$2"
  printf '%s\t%s\ttiming_only\t%s\t%s\tNA\tNA\tNA\tNA\t%s\n' \
    "${label}" \
    "${command_name}" \
    "${records}" \
    "$((records * read_len))" \
    "no normalized comparator parser configured" >> "${external_parity_tsv}"
}

stats_triplet() {
  awk -F '\t' '
    $1 == "records" { records = $2 }
    $1 == "bases" { bases = $2 }
    $1 == "checksum" { checksum = $2 }
    END { printf "%s\t%s\t%s\n", records, bases, checksum }
  '
}

microraptor_fastq_triplet() {
  target/release/microraptor stats --format fastq "$1" | stats_triplet
}

record_external_parity_triplets() {
  local label="$1"
  local tool="$2"
  local expected="$3"
  local observed="$4"
  local notes="$5"
  local expected_records expected_bases expected_checksum observed_records observed_bases observed_checksum
  IFS=$'\t' read -r expected_records expected_bases expected_checksum <<< "${expected}"
  IFS=$'\t' read -r observed_records observed_bases observed_checksum <<< "${observed}"
  local status="mismatch"
  if [[ "${expected_records}" == "${observed_records}" && "${expected_bases}" == "${observed_bases}" && "${expected_checksum}" == "${observed_checksum}" ]]; then
    status="match"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${label}" \
    "${tool}" \
    "${status}" \
    "${expected_records}" \
    "${expected_bases}" \
    "${expected_checksum}" \
    "${observed_records}" \
    "${observed_bases}" \
    "${observed_checksum}" \
    "${notes}" >> "${external_parity_tsv}"
}

record_external_parity_unknown() {
  local label="$1"
  local tool="$2"
  local notes="$3"
  printf '%s\t%s\tunknown\tNA\tNA\tNA\tNA\tNA\tNA\t%s\n' \
    "${label}" \
    "${tool}" \
    "${notes}" >> "${external_parity_tsv}"
}

record_seqkit_fastq_parity() {
  local label="$1"
  local path="$2"
  command -v seqkit >/dev/null 2>&1 || return 0
  local expected observed
  expected="$(microraptor_fastq_triplet "${path}")"
  set +e
  observed="$(seqkit seq -w 0 "${path}" 2>/dev/null | target/release/microraptor checksum --format fastq - 2>/dev/null | stats_triplet)"
  local status="$?"
  set -e
  if [[ "${status}" -ne 0 || "${observed}" == $'\t\t' ]]; then
    record_external_parity_unknown "${label}" seqkit "seqkit normalized FASTQ stream was unavailable"
    return 0
  fi
  record_external_parity_triplets "${label}" seqkit "${expected}" "${observed}" "seqkit seq -w 0 normalized FASTQ stream"
}

record_seqtk_fastq_parity() {
  local label="$1"
  local path="$2"
  command -v seqtk >/dev/null 2>&1 || return 0
  local expected observed
  expected="$(microraptor_fastq_triplet "${path}")"
  set +e
  observed="$(seqtk seq -A "${path}" 2>/dev/null | target/release/microraptor checksum --format fasta - 2>/dev/null | stats_triplet)"
  local status="$?"
  set -e
  if [[ "${status}" -ne 0 || "${observed}" == $'\t\t' ]]; then
    record_external_parity_unknown "${label}" seqtk "seqtk normalized FASTA stream was unavailable"
    return 0
  fi
  record_external_parity_triplets "${label}" seqtk "${expected}" "${observed}" "seqtk seq -A normalized FASTA sequence stream"
}

record_samtools_single_fastq_parity() {
  local label="$1"
  local path="$2"
  command -v samtools >/dev/null 2>&1 || return 0
  local expected observed
  expected="$(microraptor_fastq_triplet "${path}")"
  set +e
  observed="$(samtools import -0 "${path}" -o - -O SAM -@ "${workers}" 2>/dev/null | target/release/microraptor checksum --format sam - 2>/dev/null | stats_triplet)"
  local status="$?"
  set -e
  if [[ "${status}" -ne 0 || "${observed}" == $'\t\t' ]]; then
    record_external_parity_unknown "${label}" samtools "samtools SAM stream was unavailable"
    return 0
  fi
  record_external_parity_triplets "${label}" samtools "${expected}" "${observed}" "samtools import -O SAM sequence stream"
}

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

run_microraptor_optional() {
  local label="$1"
  local path="$2"
  [[ -f "${path}" ]] || return 0

  printf 'running microraptor %s: %s\n' "${label}" "${path}"
  set +e
  local json_output
  json_output="$(
    target/release/microraptor-bench \
      --input "${path}" \
      --iters "${iters}" \
      --workers "${workers}" \
      --json 2>&1
  )"
  local json_status="$?"
  set -e

  {
    printf '### %s\n\n' "${label}"
    printf '```text\n'
    printf '%s\n' "${json_output}"
    printf 'exit_status\t%s\n' "${json_status}"
    printf '```\n\n'
  } >> "${md}"

  if [[ "${json_status}" -eq 0 ]]; then
    printf '%s\n' "${json_output}" >> "${jsonl}"
  fi
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

run_microraptor_paired_optional() {
  local label="$1"
  local first="$2"
  local second="$3"
  [[ -f "${first}" && -f "${second}" ]] || return 0

  printf 'running microraptor %s: %s %s\n' "${label}" "${first}" "${second}"
  set +e
  local json_output
  json_output="$(
    target/release/microraptor-bench \
      --paired-inputs "${first}" "${second}" \
      --iters "${iters}" \
      --workers "${workers}" \
      --json 2>&1
  )"
  local json_status="$?"
  set -e

  {
    printf '### %s\n\n' "${label}"
    printf '```text\n'
    printf '%s\n' "${json_output}"
    printf 'exit_status\t%s\n' "${json_status}"
    printf '```\n\n'
  } >> "${md}"

  if [[ "${json_status}" -eq 0 ]]; then
    printf '%s\n' "${json_output}" >> "${jsonl}"
  fi
}

run_external() {
  local label="$1"
  local command_name="$2"
  shift 2
  local command_display="$*"

  {
    printf '### %s\n\n' "${label}"
    if command -v "${command_name}" >/dev/null 2>&1; then
      printf '```text\n'
      set +e
      local output
      output="$(/usr/bin/time -f 'elapsed_s\t%e' "$@" 2>&1)"
      local status="$?"
      set -e
      printf '%s\n' "${output}"
      local elapsed
      elapsed="$(printf '%s\n' "${output}" | awk -F '\t' '$1 == "elapsed_s" { value = $2 } END { print value }')"
      printf '%s\t%s\t%s\t%s\t%s\n' \
        "${label}" \
        "${command_name}" \
        "${status}" \
        "${elapsed:-NA}" \
        "${command_display}" >> "${external_tsv}"
      record_external_parity_timing_only "${label}" "${command_name}"
      if [[ "${status}" -ne 0 ]]; then
        printf 'exit_status\t%s\n' "${status}"
      fi
      printf '```\n\n'
    else
      printf '`%s` not installed; skipped.\n\n' "${command_name}"
      printf '%s\t%s\tskipped\tNA\t%s\n' \
        "${label}" \
        "${command_name}" \
        "${command_display}" >> "${external_tsv}"
      record_external_parity_timing_only "${label}" "${command_name}"
    fi
  } >> "${md}"
}

run_external_stdout_null() {
  local label="$1"
  local command_name="$2"
  shift 2
  local command_display="$*"

  {
    printf '### %s\n\n' "${label}"
    if command -v "${command_name}" >/dev/null 2>&1; then
      printf '```text\n'
      set +e
      local output
      output="$({ /usr/bin/time -f 'elapsed_s\t%e' "$@" >/dev/null; } 2>&1)"
      local status="$?"
      set -e
      printf '%s\n' "${output}"
      local elapsed
      elapsed="$(printf '%s\n' "${output}" | awk -F '\t' '$1 == "elapsed_s" { value = $2 } END { print value }')"
      printf '%s\t%s\t%s\t%s\t%s\n' \
        "${label}" \
        "${command_name}" \
        "${status}" \
        "${elapsed:-NA}" \
        "${command_display}" >> "${external_tsv}"
      record_external_parity_timing_only "${label}" "${command_name}"
      if [[ "${status}" -ne 0 ]]; then
        printf 'exit_status\t%s\n' "${status}"
      fi
      printf '```\n\n'
    else
      printf '`%s` not installed; skipped.\n\n' "${command_name}"
      printf '%s\t%s\tskipped\tNA\t%s\n' \
        "${label}" \
        "${command_name}" \
        "${command_display}" >> "${external_tsv}"
      record_external_parity_timing_only "${label}" "${command_name}"
    fi
  } >> "${md}"
}

safe_artifact_label() {
  printf '%s' "$1" | tr -c 'A-Za-z0-9_.-' '_'
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

if [[ "${#corpus_input_list[@]}" -gt 0 ]]; then
  for corpus_input in "${corpus_input_list[@]}"; do
    run_microraptor_optional "corpus/$(basename "${corpus_input}")" "${corpus_input}"
  done
fi

if [[ "${#corpus_paired_input_list[@]}" -gt 0 ]]; then
  for corpus_pair in "${corpus_paired_input_list[@]}"; do
    IFS=',' read -r first second label <<< "${corpus_pair}"
    label="${label:-$(basename "${first}")+$(basename "${second}")}"
    run_microraptor_paired_optional "corpus-paired/${label}" "${first}" "${second}"
  done
fi

{
  printf '## external tools\n\n'
} >> "${md}"

run_external "seqkit stats single/raw" seqkit seqkit stats "${input_dir}/single.fastq"
record_seqkit_fastq_parity "seqkit stats single/raw" "${input_dir}/single.fastq"
run_external "seqkit stats single/gzip" seqkit seqkit stats "${input_dir}/single.fastq.gz"
record_seqkit_fastq_parity "seqkit stats single/gzip" "${input_dir}/single.fastq.gz"
run_external "seqkit stats single/bgzf" seqkit seqkit stats "${input_dir}/single.fastq.bgz"
record_seqkit_fastq_parity "seqkit stats single/bgzf" "${input_dir}/single.fastq.bgz"
run_external "seqkit stats paired/r1/raw" seqkit seqkit stats "${input_dir}/r1.fastq"
record_seqkit_fastq_parity "seqkit stats paired/r1/raw" "${input_dir}/r1.fastq"
run_external "seqkit stats paired/r2/raw" seqkit seqkit stats "${input_dir}/r2.fastq"
record_seqkit_fastq_parity "seqkit stats paired/r2/raw" "${input_dir}/r2.fastq"
run_external "seqkit stats paired/r1/gzip" seqkit seqkit stats "${input_dir}/r1.fastq.gz"
record_seqkit_fastq_parity "seqkit stats paired/r1/gzip" "${input_dir}/r1.fastq.gz"
run_external "seqkit stats paired/r2/gzip" seqkit seqkit stats "${input_dir}/r2.fastq.gz"
record_seqkit_fastq_parity "seqkit stats paired/r2/gzip" "${input_dir}/r2.fastq.gz"
run_external "seqkit stats paired/r1/bgzf" seqkit seqkit stats "${input_dir}/r1.fastq.bgz"
record_seqkit_fastq_parity "seqkit stats paired/r1/bgzf" "${input_dir}/r1.fastq.bgz"
run_external "seqkit stats paired/r2/bgzf" seqkit seqkit stats "${input_dir}/r2.fastq.bgz"
record_seqkit_fastq_parity "seqkit stats paired/r2/bgzf" "${input_dir}/r2.fastq.bgz"
run_external_stdout_null "seqtk comp single/raw" seqtk seqtk comp "${input_dir}/single.fastq"
record_seqtk_fastq_parity "seqtk comp single/raw" "${input_dir}/single.fastq"
run_external_stdout_null "seqtk comp single/gzip" seqtk seqtk comp "${input_dir}/single.fastq.gz"
record_seqtk_fastq_parity "seqtk comp single/gzip" "${input_dir}/single.fastq.gz"
run_external_stdout_null "seqtk comp single/bgzf" seqtk seqtk comp "${input_dir}/single.fastq.bgz"
record_seqtk_fastq_parity "seqtk comp single/bgzf" "${input_dir}/single.fastq.bgz"
run_external "seqtk fqchk single/raw" seqtk seqtk fqchk "${input_dir}/single.fastq"
record_seqtk_fastq_parity "seqtk fqchk single/raw" "${input_dir}/single.fastq"
run_external "bgzip test single/bgzf" bgzip bgzip -t "${input_dir}/single.fastq.bgz"
run_external_stdout_null "bgzip decompress single/bgzf" bgzip bgzip -dc "${input_dir}/single.fastq.bgz"
run_external "samtools import single/raw" samtools samtools import \
  -0 "${input_dir}/single.fastq" \
  -o /dev/null \
  -O BAM \
  -@ "${workers}"
record_samtools_single_fastq_parity "samtools import single/raw" "${input_dir}/single.fastq"
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
run_external_stdout_null "fastp paired/raw" fastp fastp \
  --in1 "${input_dir}/r1.fastq" \
  --in2 "${input_dir}/r2.fastq" \
  --stdout \
  --disable_adapter_trimming \
  --disable_quality_filtering \
  --disable_length_filtering \
  --thread "${workers}" \
  --json "${result_dir}/fastp.json" \
  --html "${result_dir}/fastp.html"
run_external_stdout_null "fastp paired/gzip" fastp fastp \
  --in1 "${input_dir}/r1.fastq.gz" \
  --in2 "${input_dir}/r2.fastq.gz" \
  --stdout \
  --disable_adapter_trimming \
  --disable_quality_filtering \
  --disable_length_filtering \
  --thread "${workers}" \
  --json "${result_dir}/fastp-gzip.json" \
  --html "${result_dir}/fastp-gzip.html"
run_external_stdout_null "fastp paired/bgzf" fastp fastp \
  --in1 "${input_dir}/r1.fastq.bgz" \
  --in2 "${input_dir}/r2.fastq.bgz" \
  --stdout \
  --disable_adapter_trimming \
  --disable_quality_filtering \
  --disable_length_filtering \
  --thread "${workers}" \
  --json "${result_dir}/fastp-bgzf.json" \
  --html "${result_dir}/fastp-bgzf.html"

if [[ "${#corpus_paired_input_list[@]}" -gt 0 ]]; then
  for corpus_pair in "${corpus_paired_input_list[@]}"; do
    IFS=',' read -r first second label <<< "${corpus_pair}"
    label="${label:-$(basename "${first}")+$(basename "${second}")}"
    artifact_label="$(safe_artifact_label "${label}")"
    run_external "seqkit stats corpus-paired/${label}/r1" seqkit seqkit stats "${first}"
    run_external "seqkit stats corpus-paired/${label}/r2" seqkit seqkit stats "${second}"
    run_external "seqtk fqchk corpus-paired/${label}/r1" seqtk seqtk fqchk "${first}"
    run_external "seqtk fqchk corpus-paired/${label}/r2" seqtk seqtk fqchk "${second}"
    run_external "samtools import corpus-paired/${label}" samtools samtools import \
      -1 "${first}" \
      -2 "${second}" \
      -o /dev/null \
      -O BAM \
      -@ "${workers}"
    run_external_stdout_null "fastp corpus-paired/${label}" fastp fastp \
      --in1 "${first}" \
      --in2 "${second}" \
      --stdout \
      --disable_adapter_trimming \
      --disable_quality_filtering \
      --disable_length_filtering \
      --thread "${workers}" \
      --json "${result_dir}/fastp-${artifact_label}.json" \
      --html "${result_dir}/fastp-${artifact_label}.html"
  done
fi

if [[ "${#corpus_input_list[@]}" -gt 0 ]]; then
  for corpus_input in "${corpus_input_list[@]}"; do
    label="corpus/$(basename "${corpus_input}")"
    artifact_label="$(safe_artifact_label "${label}")"
    run_external "seqkit stats ${label}" seqkit seqkit stats "${corpus_input}"
    run_external_stdout_null "seqtk comp ${label}" seqtk seqtk comp "${corpus_input}"
    run_external "seqtk fqchk ${label}" seqtk seqtk fqchk "${corpus_input}"
    run_external "samtools import ${label}" samtools samtools import \
      -0 "${corpus_input}" \
      -o /dev/null \
      -O BAM \
      -@ "${workers}"
    run_external_stdout_null "fastp ${label}" fastp fastp \
      --in1 "${corpus_input}" \
      --stdout \
      --disable_adapter_trimming \
      --disable_quality_filtering \
      --disable_length_filtering \
      --thread "${workers}" \
      --json "${result_dir}/fastp-${artifact_label}.json" \
      --html "${result_dir}/fastp-${artifact_label}.html"
  done
fi

printf 'wrote %s\n' "${jsonl}"
printf 'wrote %s\n' "${md}"
printf 'wrote %s\n' "${metadata}"
printf 'wrote %s\n' "${external_tsv}"
