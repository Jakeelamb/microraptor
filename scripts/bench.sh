#!/usr/bin/env bash
set -euo pipefail

args=()
if [[ -n "${MICRORAPTOR_INPUT:-}" ]]; then
  args+=(--input "${MICRORAPTOR_INPUT}")
fi
if [[ -n "${MICRORAPTOR_MODE:-}" ]]; then
  args+=(--mode "${MICRORAPTOR_MODE}")
fi

cargo run --release --bin microraptor-bench -- \
  "${args[@]}" \
  --records "${MICRORAPTOR_RECORDS:-200000}" \
  --read-len "${MICRORAPTOR_READ_LEN:-150}" \
  --iters "${MICRORAPTOR_ITERS:-7}" \
  --workers "${MICRORAPTOR_WORKERS:-$(nproc)}"
