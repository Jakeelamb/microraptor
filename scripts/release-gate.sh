#!/usr/bin/env bash
set -euo pipefail

allow_dirty=0
nightly=0
bench_tier=""

usage() {
  cat <<'EOF'
usage: scripts/release-gate.sh [--allow-dirty] [--nightly] [--bench[=synthetic|fastq|fasta|all-local]]

Runs the local pre-publish gate for microraptor.

Options:
  --allow-dirty  allow a dirty git tree and pass --allow-dirty to cargo package
  --nightly      also run nightly all-feature/no-default-feature checks and fuzz build
  --bench        alias for --bench=synthetic
  --bench=fastq  regenerate synthetic FASTQ gauntlet and Rust peer snapshots
  --bench=fasta  regenerate FASTA peer and FASTA size-sweep snapshots
  --bench=all-local
                 regenerate all local synthetic FASTQ/FASTA benchmark snapshots
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
      bench_tier="synthetic"
      ;;
    --bench=*)
      bench_tier="${1#--bench=}"
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

if command -v shellcheck >/dev/null 2>&1; then
  run shellcheck --external-sources --severity=warning scripts/*.sh
else
  printf 'shellcheck not installed; CI enforces shell script linting\n' >&2
fi

if command -v cargo-deny >/dev/null 2>&1; then
  run cargo deny check
else
  printf 'cargo-deny not installed; CI enforces dependency policy\n' >&2
fi

run scripts/check-benchmark-snapshots.sh

if [[ "${nightly}" -eq 1 ]]; then
  run cargo +nightly clippy --all-targets --all-features -- -D warnings
  run cargo +nightly clippy --all-targets --no-default-features -- -D warnings
  run cargo +nightly test --all-features
  run cargo +nightly test --no-default-features
  run cargo +nightly fuzz build
fi

case "${bench_tier}" in
  "")
    ;;
  synthetic | fastq)
    run scripts/benchmark-gauntlet.sh
    run scripts/render-benchmark-report.sh
    run scripts/benchmark-rust-peers.sh
    run scripts/check-benchmark-snapshots.sh
    ;;
  fasta)
    run scripts/benchmark-fasta-peers.sh
    run scripts/benchmark-fasta-peer-size-sweep.sh
    run scripts/check-benchmark-snapshots.sh
    ;;
  all-local)
    run scripts/benchmark-gauntlet.sh
    run scripts/render-benchmark-report.sh
    run scripts/benchmark-rust-peers.sh
    run scripts/benchmark-fasta-peers.sh
    run scripts/benchmark-fasta-peer-size-sweep.sh
    run scripts/benchmark-fasta-gauntlet.sh
    run scripts/check-benchmark-snapshots.sh
    ;;
  *)
    printf 'unknown benchmark tier: %s\n' "${bench_tier}" >&2
    exit 2
    ;;
esac

run cargo package "${package_args[@]}"
