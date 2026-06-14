#!/usr/bin/env bash
set -euo pipefail

allow_dirty=0
nightly=0
bench=0

usage() {
  cat <<'EOF'
usage: scripts/release-gate.sh [--allow-dirty] [--nightly] [--bench]

Runs the local pre-publish gate for microraptor.

Options:
  --allow-dirty  allow a dirty git tree and pass --allow-dirty to cargo package
  --nightly      also run nightly all-feature/no-default-feature checks and fuzz build
  --bench        also regenerate the synthetic gauntlet and Rust peer benchmark snapshots
  -h, --help     show this help
EOF
}

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --allow-dirty)
      allow_dirty=1
      ;;
    --nightly)
      nightly=1
      ;;
    --bench)
      bench=1
      ;;
    -h | --help)
      usage
      exit 0
      ;;
    *)
      printf 'unknown argument: %s\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
  shift
done

run() {
  printf '+'
  printf ' %q' "$@"
  printf '\n'
  "$@"
}

if [[ "${allow_dirty}" -eq 0 && -n "$(git status --short)" ]]; then
  printf 'release gate requires a clean git tree; pass --allow-dirty for development verification\n' >&2
  exit 1
fi

package_args=()
if [[ "${allow_dirty}" -eq 1 ]]; then
  package_args+=(--allow-dirty)
fi

run cargo fmt --all -- --check
run env RUSTFLAGS=-D\ warnings cargo check --lib
run env RUSTDOCFLAGS=-D\ warnings cargo doc --no-deps
run cargo test --all

run test -s CITATION.cff
run test -s CHANGELOG.md
run test -s CONTRIBUTING.md
run test -s SECURITY.md
run test -s .github/pull_request_template.md
run test -s .github/ISSUE_TEMPLATE/parser_bug.yml
run test -s .github/ISSUE_TEMPLATE/benchmark_claim.yml
run test -s .github/ISSUE_TEMPLATE/feature_request.yml
run test -s docs/API_SURFACE.md
run test -s docs/FRAMEWORK.md
run test -s docs/PUBLISHING.md
run test -s docs/RELEASE_AUDIT.md

run bash -n \
  scripts/benchmark-gauntlet.sh \
  scripts/benchmark-common.sh \
  scripts/render-benchmark-report.sh \
  scripts/check-benchmark-snapshots.sh \
  scripts/benchmark-rust-peers.sh \
  scripts/benchmark-rust-peer-size-sweep.sh \
  scripts/benchmark-fasta-peers.sh \
  scripts/benchmark-fasta-peer-size-sweep.sh \
  scripts/benchmark-fasta-gauntlet.sh \
  scripts/check-replication-host.sh \
  scripts/discover-local-benchmark-corpus.sh \
  scripts/export-replication-kit.sh \
  scripts/prepare-real-benchmark-inputs.sh \
  scripts/release-gate.sh

run scripts/check-benchmark-snapshots.sh

if [[ "${nightly}" -eq 1 ]]; then
  run cargo +nightly clippy --all-targets --all-features -- -D warnings
  run cargo +nightly clippy --all-targets --no-default-features -- -D warnings
  run cargo +nightly test --all-features
  run cargo +nightly test --no-default-features
  run cargo +nightly fuzz build
fi

if [[ "${bench}" -eq 1 ]]; then
  run scripts/benchmark-gauntlet.sh
  run scripts/render-benchmark-report.sh
  run scripts/benchmark-rust-peers.sh
  run scripts/check-benchmark-snapshots.sh
fi

run cargo package "${package_args[@]}"
