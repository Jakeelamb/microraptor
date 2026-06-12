#!/usr/bin/env bash
set -euo pipefail

source_dir="${MICRORAPTOR_REAL_SOURCE_DIR:-${HOME}/Projects/Benchmarks/datasets/drosophila_melanogaster}"
out_dir="${MICRORAPTOR_REAL_INPUT_DIR:-target/bench-real-inputs/drosophila_melanogaster}"
threads="${MICRORAPTOR_REAL_BGZIP_THREADS:-$(nproc)}"

r1="${source_dir}/illumina_pe_r1.1m.fq"
r2="${source_dir}/illumina_pe_r2.1m.fq"
manifest="${out_dir}/manifest.tsv"
combined="${out_dir}/illumina_pe_r1_r2.1m-combined.fq"

display_path() {
  local path="$1"
  if [[ -n "${HOME:-}" ]]; then
    printf '%s\n' "${path}" | sed "s#^${HOME}#~#"
  else
    printf '%s\n' "${path}"
  fi
}

if [[ ! -f "${r1}" || ! -f "${r2}" ]]; then
  printf 'missing Drosophila FASTQ inputs under %s\n' "${source_dir}" >&2
  exit 1
fi
if ! command -v gzip >/dev/null 2>&1; then
  printf 'gzip is required\n' >&2
  exit 1
fi
if ! command -v bgzip >/dev/null 2>&1; then
  printf 'bgzip is required; put the benchmark conda environment on PATH\n' >&2
  exit 1
fi

mkdir -p "${out_dir}"

make_gzip() {
  local input="$1"
  local output="$2"
  if [[ ! -s "${output}" || "${input}" -nt "${output}" ]]; then
    gzip -c "${input}" > "${output}"
  fi
}

make_bgzip() {
  local input="$1"
  local output="$2"
  if [[ ! -s "${output}" || "${input}" -nt "${output}" ]]; then
    bgzip -@ "${threads}" -c "${input}" > "${output}"
  fi
}

make_gzip "${r1}" "${out_dir}/illumina_pe_r1.1m.fq.gz"
make_gzip "${r2}" "${out_dir}/illumina_pe_r2.1m.fq.gz"
make_bgzip "${r1}" "${out_dir}/illumina_pe_r1.1m.fq.bgz"
make_bgzip "${r2}" "${out_dir}/illumina_pe_r2.1m.fq.bgz"
if [[ ! -s "${combined}" || "${r1}" -nt "${combined}" || "${r2}" -nt "${combined}" ]]; then
  cat "${r1}" "${r2}" > "${combined}"
fi
make_bgzip "${combined}" "${combined}.bgz"

{
  printf 'label\tpath\tbytes\tsource\n'
  for path in \
    "${r1}" \
    "${r2}" \
    "${out_dir}/illumina_pe_r1.1m.fq.gz" \
    "${out_dir}/illumina_pe_r2.1m.fq.gz" \
    "${out_dir}/illumina_pe_r1.1m.fq.bgz" \
    "${out_dir}/illumina_pe_r2.1m.fq.bgz" \
    "${combined}" \
    "${combined}.bgz"
  do
    label="$(basename "${path}")"
    printf '%s\t%s\t%s\t%s\n' \
      "${label}" \
      "$(display_path "${path}")" \
      "$(stat -c '%s' "${path}")" \
      "$(display_path "${source_dir}")"
  done
} > "${manifest}"

printf 'wrote %s\n' "${out_dir}/illumina_pe_r1.1m.fq.gz"
printf 'wrote %s\n' "${out_dir}/illumina_pe_r2.1m.fq.gz"
printf 'wrote %s\n' "${out_dir}/illumina_pe_r1.1m.fq.bgz"
printf 'wrote %s\n' "${out_dir}/illumina_pe_r2.1m.fq.bgz"
printf 'wrote %s\n' "${combined}"
printf 'wrote %s\n' "${combined}.bgz"
printf 'wrote %s\n' "${manifest}"
