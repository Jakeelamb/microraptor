#!/usr/bin/env bash
set -euo pipefail

sizes="${MICRORAPTOR_SLAB_SIZES:-262144 1048576 4194304 8388608 16777216}"
records="${MICRORAPTOR_SLAB_RECORDS:-200000}"
read_len="${MICRORAPTOR_SLAB_READ_LEN:-150}"
iters="${MICRORAPTOR_SLAB_ITERS:-3}"
input="${MICRORAPTOR_SLAB_INPUT:-}"

extract_field() {
  local json="$1"
  local name="$2"
  local field="$3"
  sed -n "s/.*{\"name\":\"${name}\"[^}]*\"${field}\":\\([0-9]*\\).*/\\1/p" <<<"${json}" |
    head -n 1
}

best_size=""
best_ns=""

printf 'slab_size\tpack_ns\tdirect_ns\treader_ns\n'
for size in ${sizes}; do
  args=(--iters "${iters}" --mode pack --json --slab-size "${size}")
  if [[ -n "${input}" ]]; then
    args+=(--input "${input}")
  else
    args+=(--records "${records}" --read-len "${read_len}")
  fi
  json="$(cargo run --quiet --release --bin microraptor-bench -- "${args[@]}")"
  pack_ns="$(extract_field "${json}" pack-seq-qual best_ns)"
  direct_ns="$(extract_field "${json}" direct-pack-seq-qual best_ns)"
  reader_ns="$(extract_field "${json}" reader-pack-seq-qual best_ns)"
  if [[ -z "${pack_ns}" ]]; then
    pack_ns="$(extract_field "${json}" file-pack-seq-qual best_ns)"
    direct_ns="$(extract_field "${json}" file-direct-pack-seq-qual best_ns)"
    reader_ns="$(extract_field "${json}" file-reader-pack-seq-qual best_ns)"
  fi
  printf '%s\t%s\t%s\t%s\n' "${size}" "${pack_ns}" "${direct_ns}" "${reader_ns}"
  if [[ -n "${pack_ns}" && ( -z "${best_ns}" || "${pack_ns}" -lt "${best_ns}" ) ]]; then
    best_ns="${pack_ns}"
    best_size="${size}"
  fi
done

printf 'best_slab_size\t%s\n' "${best_size}"
printf 'best_pack_ns\t%s\n' "${best_ns}"
