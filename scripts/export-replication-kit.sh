#!/usr/bin/env bash
set -euo pipefail

out_dir="${1:-target/replication-kit}"
artifact_dir="${out_dir}/checked-artifacts"

case "${out_dir}" in
  "" | "/" | "." | ".." | "${HOME}" | "${HOME}/"* )
    printf 'refusing unsafe replication kit output path: %s\n' "${out_dir}" >&2
    exit 2
    ;;
  target/* )
    ;;
  * )
    if [[ "${MICRORAPTOR_ALLOW_CUSTOM_REPLICATION_OUT:-0}" != "1" ]]; then
      printf 'refusing non-target output path: %s\n' "${out_dir}" >&2
      printf 'set MICRORAPTOR_ALLOW_CUSTOM_REPLICATION_OUT=1 to override intentionally\n' >&2
      exit 2
    fi
    ;;
esac

rm -rf "${out_dir}"
mkdir -p "${artifact_dir}"

copy_if_present() {
  local path="$1"
  if [[ -e "${path}" ]]; then
    mkdir -p "${artifact_dir}/$(dirname "${path}")"
    cp -R "${path}" "${artifact_dir}/${path}"
  fi
}

for path in \
  README.md \
  BENCHMARKING.md \
  CHANGELOG.md \
  CITATION.cff \
  docs/API_SURFACE.md \
  docs/FRAMEWORK.md \
  docs/PUBLISHING.md \
  docs/RELEASE_AUDIT.md \
  docs/REPLICATION.md \
  docs/benchmarks/CORPUS.md \
  scripts/benchmark-common.sh \
  scripts/benchmark-gauntlet.sh \
  scripts/render-benchmark-report.sh \
  scripts/check-benchmark-snapshots.sh \
  scripts/benchmark-rust-peers.sh \
  scripts/benchmark-rust-peer-size-sweep.sh \
  scripts/benchmark-fasta-peers.sh \
  scripts/benchmark-fasta-peer-size-sweep.sh \
  scripts/benchmark-fasta-gauntlet.sh \
  scripts/check-replication-host.sh \
  scripts/discover-local-benchmark-corpus.sh \
  scripts/prepare-real-benchmark-inputs.sh \
  scripts/release-gate.sh; do
  copy_if_present "${path}"
done

for snapshot in docs/benchmarks/*; do
  [[ -d "${snapshot}" ]] || continue
  copy_if_present "${snapshot}"
done

{
  printf '# Microraptor Replication Kit\n\n'
  printf -- '- generated_at_utc: %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  printf -- '- git_commit: %s\n' "$(git rev-parse --short HEAD 2>/dev/null || printf unknown)"
  printf -- '- git_dirty: %s\n' "$(if [[ -n "$(git status --short 2>/dev/null)" ]]; then printf true; else printf false; fi)"
  printf -- '- rustc: %s\n' "$(rustc --version 2>/dev/null || printf unavailable)"
  printf -- '- cargo: %s\n' "$(cargo --version 2>/dev/null || printf unavailable)"
  printf -- '- kernel: %s\n' "$(uname -srmo 2>/dev/null || printf unavailable)"
  printf '\n'
  printf 'This kit contains checked benchmark summaries, raw JSONL/TSV artifacts,\n'
  printf 'figures, release documentation, and scripts needed to reproduce or compare\n'
  printf 'microraptor benchmark evidence on another machine. It does not contain FASTQ\n'
  printf 'datasets; see `checked-artifacts/docs/benchmarks/CORPUS.md` and\n'
  printf '`checked-artifacts/BENCHMARKING.md` for dataset acquisition and local corpus\n'
  printf 'discovery rules.\n\n'
  printf '## Minimum Verification\n\n'
  printf '```bash\n'
  printf 'scripts/release-gate.sh\n'
  printf 'scripts/check-replication-host.sh --strict\n'
  printf 'scripts/check-benchmark-snapshots.sh\n'
  printf '```\n\n'
  printf '## Independent-Machine Benchmark Commands\n\n'
  printf '```bash\n'
  printf 'scripts/benchmark-gauntlet.sh\n'
  printf 'scripts/render-benchmark-report.sh target/bench-results/microraptor-gauntlet.jsonl docs/benchmarks/latest\n'
  printf 'MICRORAPTOR_RUST_PEER_ITERS=3 scripts/benchmark-rust-peers.sh\n'
  printf 'PATH=~/miniconda3/envs/bench/bin:$PATH scripts/discover-local-benchmark-corpus.sh\n'
  printf 'source target/bench-corpus/recommended-gauntlet.env\n'
  printf 'PATH=~/miniconda3/envs/bench/bin:$PATH MICRORAPTOR_GAUNTLET_RESULT_DIR=target/bench-results-real scripts/benchmark-gauntlet.sh\n'
  printf 'scripts/render-benchmark-report.sh target/bench-results-real/microraptor-gauntlet.jsonl docs/benchmarks/real-snapshot\n'
  printf 'scripts/check-benchmark-snapshots.sh\n'
  printf '```\n\n'
  printf '## Claim Boundary\n\n'
  printf 'Treat checked artifacts as machine-local evidence. A broad performance claim\n'
  printf 'requires regenerated artifacts from a clean release commit plus at least one\n'
  printf 'independent-machine replication run with hardware, toolchain, comparator\n'
  printf 'versions, raw JSONL/TSV, summaries, and figures.\n'
} > "${out_dir}/README.md"

{
  printf '# Replication Artifact Checksums\n\n'
  printf 'Generated from checked artifacts only. Paths are relative to the replication\n'
  printf 'kit root.\n\n'
  printf '```text\n'
  (
    cd "${out_dir}"
    find checked-artifacts -type f -print0 | sort -z | xargs -0 sha256sum
  )
  printf '```\n'
} > "${out_dir}/SHA256SUMS.md"

printf 'wrote %s\n' "${out_dir}/README.md"
printf 'wrote %s\n' "${out_dir}/SHA256SUMS.md"
