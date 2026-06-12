#!/usr/bin/env bash
set -euo pipefail

strict=0

usage() {
  cat <<'EOF'
usage: scripts/check-replication-host.sh [--strict]

Checks whether the local host has the tools needed to regenerate microraptor
release and benchmark evidence. The script does not install anything.

Options:
  --strict       require benchmark comparators and nightly/fuzz tooling too
  -h, --help     show this help
EOF
}

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --strict)
      strict=1
      ;;
    -h | --help)
      usage
      exit 0
      ;;
    -*)
      printf 'unknown argument: %s\n' "$1" >&2
      usage >&2
      exit 2
      ;;
    *)
      printf 'unexpected argument: %s\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
  shift
done

host_label="$(hostname 2>/dev/null || printf unknown)"
missing_core=0
missing_strict=0

one_line() {
  tr '\t\n' '  ' | sed 's/[[:space:]][[:space:]]*/ /g; s/[[:space:]]$//'
}

print_row() {
  local status="$1"
  local requirement="$2"
  local tool="$3"
  local detail="$4"
  printf '%s\t%s\t%s\t%s\t%s\n' "${status}" "${host_label}" "${requirement}" "${tool}" "${detail}"
}

version_for() {
  local tool="$1"
  case "${tool}" in
    bash) { bash --version 2>&1 || true; } | head -n 1 | one_line ;;
    git) { git --version 2>&1 || true; } | head -n 1 | one_line ;;
    jq) { jq --version 2>&1 || true; } | head -n 1 | one_line ;;
    cargo) { cargo --version 2>&1 || true; } | head -n 1 | one_line ;;
    rustc) { rustc --version 2>&1 || true; } | head -n 1 | one_line ;;
    seqkit) { seqkit version 2>&1 || true; } | head -n 1 | one_line ;;
    seqtk) { seqtk 2>&1 || true; } | head -n 1 | one_line ;;
    samtools) { samtools --version 2>&1 || true; } | head -n 1 | one_line ;;
    bgzip) { bgzip --version 2>&1 || true; } | head -n 1 | one_line ;;
    fastp) { fastp --version 2>&1 || true; } | head -n 1 | one_line ;;
    *) { "${tool}" --version 2>&1 || true; } | head -n 1 | one_line ;;
  esac
}

check_command() {
  local requirement="$1"
  local tool="$2"
  local required="$3"
  local path=''
  local detail=''

  path="$(command -v "${tool}" 2>/dev/null || true)"
  if [[ -n "${path}" ]]; then
    detail="$(version_for "${tool}")"
    print_row ok "${requirement}" "${tool}" "${detail} (${path})"
  else
    print_row missing "${requirement}" "${tool}" unavailable
    if [[ "${required}" == core ]]; then
      missing_core=$((missing_core + 1))
    elif [[ "${required}" == strict && "${strict}" -eq 1 ]]; then
      missing_strict=$((missing_strict + 1))
    fi
  fi
}

check_cargo_subcommand() {
  local requirement="$1"
  local label="$2"
  local command_args=("${@:3}")
  local detail=''

  if detail="$("${command_args[@]}" 2>&1 | head -n 1 | one_line)"; then
    print_row ok "${requirement}" "${label}" "${detail}"
  else
    print_row missing "${requirement}" "${label}" "${detail:-unavailable}"
    if [[ "${strict}" -eq 1 ]]; then
      missing_strict=$((missing_strict + 1))
    fi
  fi
}

printf 'status\thost\trequirement\ttool\tdetail\n'
print_row info platform hostname "${host_label}"
print_row info platform kernel "$(uname -srmo 2>/dev/null | one_line || printf unavailable)"
print_row info platform filesystem "$(df -h . 2>/dev/null | tail -n 1 | one_line || printf unavailable)"
print_row info platform cpu "$(awk -F ': ' '/model name/ { print $2; exit }' /proc/cpuinfo 2>/dev/null | one_line || printf unavailable)"
print_row info platform memory "$(awk '/MemTotal/ { printf "%.1f GiB", $2 / 1024 / 1024 }' /proc/meminfo 2>/dev/null || printf unavailable)"

check_command core bash core
check_command core git core
check_command core jq core
check_command core cargo core
check_command core rustc core

check_cargo_subcommand nightly cargo+nightly cargo +nightly --version
check_cargo_subcommand fuzz cargo-fuzz cargo fuzz --version

check_command comparator seqkit strict
check_command comparator seqtk strict
check_command comparator samtools strict
check_command comparator bgzip strict
check_command comparator fastp strict

if [[ "${missing_core}" -gt 0 ]]; then
  printf 'missing %s core replication prerequisite(s)\n' "${missing_core}" >&2
  exit 1
fi

if [[ "${strict}" -eq 1 && "${missing_strict}" -gt 0 ]]; then
  printf 'missing %s strict benchmark prerequisite(s)\n' "${missing_strict}" >&2
  exit 1
fi
