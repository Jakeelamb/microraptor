#!/usr/bin/env bash
set -euo pipefail

scripts/profile-hotpath.sh >/dev/null

profile_dir="${MICRORAPTOR_PROFILE_DIR:-target/profiles}"
jsonl="${profile_dir}/microraptor-hotpath.jsonl"
perf_stat="${profile_dir}/microraptor-pack.perf-stat.txt"

if [[ ! -f "${perf_stat}" ]]; then
  printf 'perf stat output missing; run on a host with perf permissions\n' >&2
  exit 2
fi

bases="$(
  sed -n 's/.*"name":"file-pack-seq-qual"[^}]*"bases":\([0-9]*\).*/\1/p' "${jsonl}" |
    tail -n 1
)"
instructions="$(
  awk '/instructions:u/ { gsub(/[^0-9]/, "", $1); print $1; exit }' "${perf_stat}"
)"

if [[ -z "${bases}" || -z "${instructions}" ]]; then
  printf 'could not extract bases/instructions\n' >&2
  exit 2
fi

scaled="$(( instructions * 100 / bases ))"
whole="$(( scaled / 100 ))"
frac="$(( scaled % 100 ))"
printf 'pack_instructions\t%s\n' "${instructions}"
printf 'pack_bases\t%s\n' "${bases}"
printf 'pack_instructions_per_base\t%d.%02d\n' "${whole}" "${frac}"

if [[ -n "${MICRORAPTOR_MAX_PACK_INSTRUCTIONS_PER_BASE:-}" ]]; then
  max_scaled="$(awk -v max="${MICRORAPTOR_MAX_PACK_INSTRUCTIONS_PER_BASE}" 'BEGIN { printf "%d", max * 100 }')"
  if (( scaled > max_scaled )); then
    printf 'instruction regression: %d.%02d > %s\n' \
      "${whole}" "${frac}" "${MICRORAPTOR_MAX_PACK_INSTRUCTIONS_PER_BASE}" >&2
    exit 1
  fi
fi
