#!/usr/bin/env bash
set -euo pipefail

mkdir -p target/profiles

cargo build --release --bin microraptor-bench
perf stat -d -r "${MICRORAPTOR_PERF_REPEATS:-3}" \
  target/release/microraptor-bench \
  --records "${MICRORAPTOR_RECORDS:-500000}" \
  --read-len "${MICRORAPTOR_READ_LEN:-150}" \
  --iters "${MICRORAPTOR_ITERS:-3}" \
  --workers "${MICRORAPTOR_WORKERS:-$(nproc)}"

perf record -F 999 -g -o target/profiles/microraptor-bench.perf.data -- \
  target/release/microraptor-bench \
  --records "${MICRORAPTOR_RECORDS:-500000}" \
  --read-len "${MICRORAPTOR_READ_LEN:-150}" \
  --iters 1 \
  --workers "${MICRORAPTOR_WORKERS:-$(nproc)}"

perf report -i target/profiles/microraptor-bench.perf.data --stdio \
  > target/profiles/microraptor-bench.perf.txt

echo "wrote target/profiles/microraptor-bench.perf.data"
echo "wrote target/profiles/microraptor-bench.perf.txt"
