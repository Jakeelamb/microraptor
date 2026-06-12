#!/usr/bin/env bash
set -euo pipefail

cargo run --release --bin microraptor-bench -- \
  --records "${MICRORAPTOR_RECORDS:-200000}" \
  --read-len "${MICRORAPTOR_READ_LEN:-150}" \
  --iters "${MICRORAPTOR_ITERS:-7}" \
  --workers "${MICRORAPTOR_WORKERS:-$(nproc)}"
