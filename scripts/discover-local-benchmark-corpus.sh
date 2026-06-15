#!/usr/bin/env bash
set -euo pipefail

bench_dir="${MICRORAPTOR_BENCHMARKS_DIR:-${HOME}/Projects/Benchmarks}"
out_dir="${MICRORAPTOR_CORPUS_OUT_DIR:-target/bench-corpus}"
dros_dir="${bench_dir}/datasets/drosophila_melanogaster"
prepared="${bench_dir}/manifests/drosophila_prepared.tsv"
local_manifest="${bench_dir}/manifests/local_datasets.tsv"
out_tsv="${out_dir}/local-corpus.tsv"
env_file="${out_dir}/recommended-gauntlet.env"
independent_env_file="${out_dir}/independent-gauntlet.env"
larger_env_file="${out_dir}/larger-gauntlet.env"

mkdir -p "${out_dir}"

display_path() {
  local path="$1"
  if [[ -n "${HOME:-}" ]]; then
    printf '%s\n' "${path}" | sed "s#^${HOME}#~#"
  else
    printf '%s\n' "${path}"
  fi
}

sha256_file() {
  if [[ -n "$1" && -f "$1" ]]; then
    sha256sum "$1" | awk '{ print $1 }'
  else
    printf '\n'
  fi
}

prepared_field() {
  local prepared_id="$1"
  local field="$2"
  awk -F '\t' -v id="${prepared_id}" -v field="${field}" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == field) {
          col = i
        }
      }
      next
    }
    $1 == id && col {
      print $col
      found = 1
      exit
    }
    END {
      if (!found) {
        exit 1
      }
    }
  ' "${prepared}" 2>/dev/null || true
}

scale_to_records() {
  local scale="$1"
  case "${scale}" in
    1m) printf '1000000' ;;
    5m) printf '5000000' ;;
    25m) printf '25000000' ;;
    50m) printf '50000000' ;;
    *) printf 'unknown' ;;
  esac
}

scale_rank() {
  case "$1" in
    1m) printf '1' ;;
    5m) printf '2' ;;
    25m) printf '3' ;;
    50m) printf '4' ;;
    *) printf '99' ;;
  esac
}

read_manifest_value() {
  local id="$1"
  local field="$2"
  if [[ -f "${prepared}" ]]; then
    prepared_field "${id}" "${field}"
  fi
}

join_lines() {
  local IFS=$'\n'
  printf '%s' "$*"
}

write_exported_list() {
  local name="$1"
  shift
  local value
  value="$(join_lines "$@")"
  if [[ "$#" -gt 0 ]]; then
    value+=$'\n'
  fi
  printf 'export %s=%q\n' "${name}" "${value}"
}

printf 'label\tread_type\tlayout\trole\trecords\tbases\tpath_a\tsha256_a\tpath_b\tsha256_b\tsource_manifest\tnotes\n' > "${out_tsv}"

if [[ ! -d "${bench_dir}" ]]; then
  printf 'missing benchmark workspace: %s\n' "${bench_dir}" >&2
  exit 1
fi

if [[ -d "${dros_dir}" ]]; then
  for id in pacbio_clr_subreads_50k ont_50k; do
    rel_path="$(read_manifest_value "${id}" path)"
    records="$(read_manifest_value "${id}" records)"
    bases="$(read_manifest_value "${id}" bases)"
    if [[ -n "${rel_path}" && -f "${bench_dir}/${rel_path}" ]]; then
      case "${id}" in
        pacbio_clr_subreads_50k) read_type='PacBio CLR/subreads' ;;
        ont_50k) read_type='Oxford Nanopore' ;;
        *) read_type='single-end FASTQ' ;;
      esac
      printf '%s\t%s\tsingle-end\tread-type-coverage\t%s\t%s\t%s\t%s\t\t\t%s\t%s\n' \
        "drosophila_${id}" \
        "${read_type}" \
        "${records:-unknown}" \
        "${bases:-unknown}" \
        "$(display_path "${bench_dir}/${rel_path}")" \
        "$(sha256_file "${bench_dir}/${rel_path}")" \
        "$(display_path "${prepared}")" \
        'real Drosophila long-read FASTQ; parser evidence only, not assembly quality evidence' >> "${out_tsv}"
    fi
  done

  for scale in 1m 5m 25m 50m; do
    r1="${dros_dir}/illumina_pe_r1.${scale}.fq"
    [[ -f "${r1}" ]] || continue
    r2="${dros_dir}/illumina_pe_r2.${scale}.fq"
    [[ -f "${r2}" ]] || continue

    id1="illumina_pe_r1_${scale}"
    id2="illumina_pe_r2_${scale}"
    r1_records="$(read_manifest_value "${id1}" records)"
    r2_records="$(read_manifest_value "${id2}" records)"
    r1_bases="$(read_manifest_value "${id1}" bases)"
    r2_bases="$(read_manifest_value "${id2}" bases)"

    if [[ -n "${r1_records}" && -n "${r2_records}" ]]; then
      records="$((r1_records + r2_records))"
    else
      records="$(scale_to_records "${scale}")"
      if [[ "${records}" != unknown ]]; then
        records="$((records * 2))"
      fi
    fi

    if [[ -n "${r1_bases}" && -n "${r2_bases}" ]]; then
      bases="$((r1_bases + r2_bases))"
    else
      bases='unknown'
    fi

    rank="$(scale_rank "${scale}")"
    if [[ "${rank}" -le 1 ]]; then
      role='recommended-release-corpus'
      notes='bounded paired Illumina rung for routine local release evidence'
    elif [[ "${rank}" -eq 2 ]]; then
      role='larger-release-corpus'
      notes='larger paired Illumina rung; run explicitly because external workflow comparators can dominate wall time'
    else
      role='stress-corpus'
      notes='large paired Illumina rung; run deliberately because wall time and disk pressure are higher'
    fi

    printf '%s\tIllumina paired-end\tpaired\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "drosophila_illumina_${scale}_paired_raw" \
      "${role}" \
      "${records}" \
      "${bases}" \
      "$(display_path "${r1}")" \
      "$(sha256_file "${r1}")" \
      "$(display_path "${r2}")" \
      "$(sha256_file "${r2}")" \
      "$(display_path "${prepared}")" \
      "${notes}" >> "${out_tsv}"
  done
fi

if [[ -f "${local_manifest}" ]]; then
  while IFS=$'\t' read -r dataset_id _organism read_type scale local_path _reference notes; do
    [[ "${dataset_id}" != "dataset_id" ]] || continue
    IFS=';' read -r first_rel second_rel <<< "${local_path}"
    [[ -n "${first_rel}" && -n "${second_rel}" ]] || continue
    first="${bench_dir}/${first_rel}"
    second="${bench_dir}/${second_rel}"
    [[ "${first}" != "${bench_dir}/datasets/drosophila_melanogaster/"* ]] || continue
    [[ -f "${first}" && -f "${second}" ]] || continue

    records="unknown"
    case "${scale}" in
      *"100000 read pairs"*) records=200000 ;;
      *"10000 read pairs"*) records=20000 ;;
      *"1000 read pairs"*) records=2000 ;;
    esac
    role="smoke-corpus"
    case "${dataset_id}" in
      ecoli_mg1655_srr001666_100k_pairs|yeast_btt_err1308583_10k_pairs)
        role="independent-release-corpus"
        ;;
    esac
    manifest="${bench_dir}/manifests/local_datasets.tsv"
    printf '%s\t%s\tpaired\t%s\t%s\tunknown\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${dataset_id}" \
      "${read_type}" \
      "${role}" \
      "${records}" \
      "$(display_path "${first}")" \
      "$(sha256_file "${first}")" \
      "$(display_path "${second}")" \
      "$(sha256_file "${second}")" \
      "$(display_path "${manifest}")" \
      "${notes}" >> "${out_tsv}"
  done < "${local_manifest}"
fi

recommended_inputs=()
recommended_pairs=()
independent_pairs=()
larger_pairs=()
if [[ -f "${dros_dir}/pacbio_clr_subreads.50k.fq" ]]; then
  recommended_inputs+=("${dros_dir}/pacbio_clr_subreads.50k.fq")
fi
if [[ -f "${dros_dir}/ont.50k.fq" ]]; then
  recommended_inputs+=("${dros_dir}/ont.50k.fq")
fi
scale="1m"
r1="${dros_dir}/illumina_pe_r1.${scale}.fq"
r2="${dros_dir}/illumina_pe_r2.${scale}.fq"
if [[ -f "${r1}" && -f "${r2}" ]]; then
  recommended_pairs+=("${r1},${r2},drosophila_illumina_${scale}_raw")
fi
scale="5m"
r1="${dros_dir}/illumina_pe_r1.${scale}.fq"
r2="${dros_dir}/illumina_pe_r2.${scale}.fq"
if [[ -f "${r1}" && -f "${r2}" ]]; then
  larger_pairs+=("${r1},${r2},drosophila_illumina_${scale}_raw")
fi
if [[ -f "${local_manifest}" ]]; then
  while IFS=$'\t' read -r dataset_id _organism _read_layout _scale local_path _reference _notes; do
    case "${dataset_id}" in
      ecoli_mg1655_srr001666_100k_pairs|yeast_btt_err1308583_10k_pairs)
        IFS=';' read -r first second <<< "${local_path}"
        first="${bench_dir}/${first}"
        second="${bench_dir}/${second}"
        if [[ -f "${first}" && -f "${second}" ]]; then
          independent_pairs+=("${first},${second},${dataset_id}")
        fi
        ;;
    esac
  done < <(tail -n +2 "${local_manifest}")
fi

{
  printf '# shellcheck shell=bash\n'
  printf '# Generated by scripts/discover-local-benchmark-corpus.sh\n'
  printf '# Corpus variables are newline-delimited so paths may contain spaces.\n'
  write_exported_list MICRORAPTOR_GAUNTLET_CORPUS_INPUTS "${recommended_inputs[@]}"
  write_exported_list MICRORAPTOR_GAUNTLET_CORPUS_PAIRED_INPUTS "${recommended_pairs[@]}"
} > "${env_file}"

{
  printf '# shellcheck shell=bash\n'
  printf '# Generated by scripts/discover-local-benchmark-corpus.sh\n'
  printf '# Independent non-Drosophila local corpus rows for replication evidence.\n'
  printf '# Corpus variables are newline-delimited so paths may contain spaces.\n'
  write_exported_list MICRORAPTOR_GAUNTLET_CORPUS_INPUTS
  write_exported_list MICRORAPTOR_GAUNTLET_CORPUS_PAIRED_INPUTS "${independent_pairs[@]}"
} > "${independent_env_file}"

{
  printf '# shellcheck shell=bash\n'
  printf '# Generated by scripts/discover-local-benchmark-corpus.sh\n'
  printf '# Source after recommended-gauntlet.env to add larger raw paired rows.\n'
  printf '# Corpus variables are newline-delimited so paths may contain spaces.\n'
  if [[ "${#larger_pairs[@]}" -gt 0 ]]; then
    printf 'if [[ -n "${MICRORAPTOR_GAUNTLET_CORPUS_PAIRED_INPUTS:-}" ]]; then\n'
    printf "  MICRORAPTOR_GAUNTLET_CORPUS_PAIRED_INPUTS+=\$'\\n'\n"
    printf 'fi\n'
    printf 'MICRORAPTOR_GAUNTLET_CORPUS_PAIRED_INPUTS+=%q\n' "$(join_lines "${larger_pairs[@]}")"$'\n'
    printf 'export MICRORAPTOR_GAUNTLET_CORPUS_PAIRED_INPUTS\n'
  else
    printf 'export MICRORAPTOR_GAUNTLET_CORPUS_PAIRED_INPUTS="${MICRORAPTOR_GAUNTLET_CORPUS_PAIRED_INPUTS:-}"\n'
  fi
} > "${larger_env_file}"

printf 'wrote %s\n' "${out_tsv}"
printf 'wrote %s\n' "${env_file}"
printf 'wrote %s\n' "${independent_env_file}"
printf 'wrote %s\n' "${larger_env_file}"
