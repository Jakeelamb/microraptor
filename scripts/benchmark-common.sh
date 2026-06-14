#!/usr/bin/env bash

microraptor_set_thread_cap() {
  local bench_threads="$1"
  export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-${bench_threads}}"
  export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-${bench_threads}}"
  export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${bench_threads}}"
  export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-${bench_threads}}"
  export MKL_NUM_THREADS="${MKL_NUM_THREADS:-${bench_threads}}"
}

microraptor_display_path() {
  if [[ -n "${HOME:-}" ]]; then
    awk -v home="${HOME}" 'BEGIN { value = ARGV[1]; ARGV[1] = ""; gsub(home, "~", value); print value }' "$1"
  else
    printf '%s\n' "$1"
  fi
}

microraptor_write_sanitized_file() {
  local src="$1"
  local dst="$2"
  local tmp
  tmp="$(mktemp "${dst}.tmp.XXXXXX")"
  if [[ -n "${HOME:-}" ]]; then
    awk -v home="${HOME}" '{ gsub(home, "~"); print }' "${src}" > "${tmp}"
  else
    cp "${src}" "${tmp}"
  fi
  mv "${tmp}" "${dst}"
}

microraptor_feature_list() {
  local features="$1"
  printf '%s\n' "${features}" | awk -F ',' '{
    for (i = 1; i <= NF; i++) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $i)
      if ($i != "") {
        printf "%s\"%s\"", sep, $i
        sep = ", "
      }
    }
  }'
}

microraptor_dependency_spec() {
  local features="$1"
  local feature_list
  if [[ -n "${features}" ]]; then
    feature_list="$(microraptor_feature_list "${features}")"
    printf 'microraptor = { path = "../..", features = [%s] }\n' "${feature_list}"
  else
    printf 'microraptor = { path = "../.." }\n'
  fi
}
