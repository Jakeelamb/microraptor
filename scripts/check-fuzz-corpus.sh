#!/usr/bin/env bash
set -euo pipefail

if ! command -v cargo-fuzz >/dev/null 2>&1 && ! cargo +nightly fuzz --help >/dev/null 2>&1; then
  printf 'cargo-fuzz is required for fuzz corpus smoke tests\n' >&2
  exit 1
fi

shopt -s nullglob
checked=0
for corpus_dir in fuzz/corpus/*; do
  [[ -d "${corpus_dir}" ]] || continue
  target="$(basename "${corpus_dir}")"
  ASAN_OPTIONS="${ASAN_OPTIONS:-detect_leaks=0}" \
    LSAN_OPTIONS="${LSAN_OPTIONS:-detect_leaks=0}" \
    cargo +nightly fuzz run "${target}" "${corpus_dir}" -- -runs=0
  checked=$((checked + 1))
done

if [[ "${checked}" -eq 0 ]]; then
  printf 'no fuzz corpus directories found\n' >&2
  exit 1
fi

printf 'checked %s fuzz corpus target(s)\n' "${checked}"
