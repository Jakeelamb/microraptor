#!/usr/bin/env bash
set -euo pipefail

records="${MICRORAPTOR_PROFILE_RECORDS:-500000}"
read_len="${MICRORAPTOR_PROFILE_READ_LEN:-150}"
iters="${MICRORAPTOR_PROFILE_ITERS:-3}"
workers="${MICRORAPTOR_WORKERS:-$(nproc)}"
input_dir="${MICRORAPTOR_PROFILE_INPUT_DIR:-target/bench-inputs}"
profile_dir="${MICRORAPTOR_PROFILE_DIR:-target/profiles}"

mkdir -p "${input_dir}" "${profile_dir}"

cargo build --release --bin microraptor-bench --bin microraptor-fixture
target/release/microraptor-fixture \
  --out-dir "${input_dir}" \
  --records "${records}" \
  --read-len "${read_len}"

input="${MICRORAPTOR_PROFILE_INPUT:-${input_dir}/single.fastq}"
jsonl="${profile_dir}/microraptor-hotpath.jsonl"
: > "${jsonl}"

for mode in parse pack; do
  printf 'benchmarking mode=%s input=%s\n' "${mode}" "${input}"
  target/release/microraptor-bench \
    --input "${input}" \
    --mode "${mode}" \
    --iters "${iters}" \
    --workers "${workers}" \
    --json >> "${jsonl}"

  if command -v perf >/dev/null 2>&1; then
    perf stat -d \
      -o "${profile_dir}/microraptor-${mode}.perf-stat.txt" \
      target/release/microraptor-bench \
        --input "${input}" \
        --mode "${mode}" \
        --iters "${iters}" \
        --workers "${workers}" >/dev/null || true

    perf record -F 999 -g \
      -o "${profile_dir}/microraptor-${mode}.perf.data" -- \
      target/release/microraptor-bench \
        --input "${input}" \
        --mode "${mode}" \
        --iters 1 \
        --workers "${workers}" >/dev/null 2>&1 || true

    if [[ -f "${profile_dir}/microraptor-${mode}.perf.data" ]]; then
      perf report \
        -i "${profile_dir}/microraptor-${mode}.perf.data" \
        --stdio > "${profile_dir}/microraptor-${mode}.perf-report.txt" || true
    fi
  fi
done

printf 'wrote %s\n' "${jsonl}"
printf 'wrote %s/microraptor-parse.perf-stat.txt if perf was permitted\n' "${profile_dir}"
printf 'wrote %s/microraptor-pack.perf-stat.txt if perf was permitted\n' "${profile_dir}"
