#!/usr/bin/env bash
set -euo pipefail

records="${MICRORAPTOR_PACK_REGRESSION_RECORDS:-500000}"
read_len="${MICRORAPTOR_PACK_REGRESSION_READ_LEN:-150}"
iters="${MICRORAPTOR_PACK_REGRESSION_ITERS:-7}"
tolerance_pct="${MICRORAPTOR_PACK_REGRESSION_TOLERANCE_PCT:-15}"

json="$(
  cargo run --release --bin microraptor-bench -- \
    --records "${records}" \
    --read-len "${read_len}" \
    --iters "${iters}" \
    --mode pack \
    --json
)"

extract_field() {
  local name="$1"
  local field="$2"
  printf '%s\n' "${json}" |
    sed -n "s/.*{\"name\":\"${name}\"[^}]*\"${field}\":\\([0-9]*\\).*/\\1/p"
}

fast_ns="$(extract_field pack-seq-qual best_ns)"
reference_ns="$(extract_field reader-pack-seq-qual best_ns)"
fast_checksum="$(extract_field pack-seq-qual checksum)"
reference_checksum="$(extract_field reader-pack-seq-qual checksum)"

if [[ -z "${fast_ns}" || -z "${reference_ns}" ]]; then
  printf 'missing pack benchmark rows\n%s\n' "${json}" >&2
  exit 2
fi

if [[ "${fast_checksum}" != "${reference_checksum}" ]]; then
  printf 'checksum mismatch: pack=%s reader=%s\n' "${fast_checksum}" "${reference_checksum}" >&2
  exit 1
fi

limit_ns="$(( reference_ns + (reference_ns * tolerance_pct / 100) ))"
printf 'pack-seq-qual_ns\t%s\n' "${fast_ns}"
printf 'reader-pack-seq-qual_ns\t%s\n' "${reference_ns}"
printf 'tolerance_pct\t%s\n' "${tolerance_pct}"

if (( fast_ns > limit_ns )); then
  printf 'pack regression: %s ns > %s ns\n' "${fast_ns}" "${limit_ns}" >&2
  exit 1
fi
