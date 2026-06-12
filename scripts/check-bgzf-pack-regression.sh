#!/usr/bin/env bash
set -euo pipefail

records="${MICRORAPTOR_BGZF_PACK_REGRESSION_RECORDS:-100000}"
read_len="${MICRORAPTOR_BGZF_PACK_REGRESSION_READ_LEN:-150}"
iters="${MICRORAPTOR_BGZF_PACK_REGRESSION_ITERS:-7}"
large_records="${MICRORAPTOR_BGZF_PACK_REGRESSION_LARGE_RECORDS:-300000}"
large_iters="${MICRORAPTOR_BGZF_PACK_REGRESSION_LARGE_ITERS:-3}"
workers="${MICRORAPTOR_WORKERS:-$(nproc)}"
tolerance_pct="${MICRORAPTOR_BGZF_PACK_REGRESSION_TOLERANCE_PCT:-15}"
parallel_min_bytes="${MICRORAPTOR_BGZF_PARALLEL_MIN_BYTES:-33554432}"
input_dir="${MICRORAPTOR_BGZF_PACK_REGRESSION_INPUT_DIR:-target/bench-inputs}"
features="${MICRORAPTOR_BGZF_PACK_FEATURES:-libdeflate}"
timing_args=()
if [[ "${CI:-}" == "true" && "${MICRORAPTOR_ENFORCE_TIMING:-0}" != "1" ]]; then
  timing_args+=(--skip-timing-checks)
fi

cargo_run=(cargo run --quiet --release)
if [[ -n "${features}" ]]; then
  cargo_run+=(--features "${features}")
fi

check_case() {
  local label="$1"
  local dir="$2"
  local case_records="$3"
  local case_iters="$4"
  local pattern="$5"
  local min_input_bytes="$6"

  "${cargo_run[@]}" --bin microraptor-fixture -- \
    --out-dir "${dir}" \
    --records "${case_records}" \
    --read-len "${read_len}" \
    --pattern "${pattern}" >/dev/null

  "${cargo_run[@]}" --bin microraptor-bench -- \
    --input "${dir}/single.fastq.bgz" \
    --iters "${case_iters}" \
    --workers "${workers}" \
    --mode pack \
    --check-bgzf-pack-regression \
    --check-label "${label}" \
    --min-input-bytes "${min_input_bytes}" \
    --tolerance-pct "${tolerance_pct}" \
    "${timing_args[@]}"
}

check_case small "${input_dir}/bgzf-small" "${records}" "${iters}" cyclic 0
if (( large_records > 0 )); then
  check_case large "${input_dir}/bgzf-large" "${large_records}" "${large_iters}" entropy "${parallel_min_bytes}"
fi
printf 'tolerance_pct\t%s\n' "${tolerance_pct}"
