#!/usr/bin/env bash
set -euo pipefail

out_dir="${MICRORAPTOR_ASM_DIR:-target/asm}"
mkdir -p "${out_dir}"

cargo rustc --release --all-features --lib -- --emit=asm

asm_file="$(
  find target/release/deps -maxdepth 1 -type f -name 'microraptor-*.s' \
    -printf '%T@ %p\n' |
    sort -nr |
    awk 'NR == 1 { print $2 }'
)"

if [[ -z "${asm_file}" ]]; then
  printf 'no assembly file emitted\n' >&2
  exit 1
fi

dest="${out_dir}/microraptor-pack.s"
cp "${asm_file}" "${dest}"

printf 'asm\t%s\n' "${dest}"
printf 'kernel_symbols\t'
grep -E 'pack_trusted_fastq|pack_bases_and_qualities|summarize_qualities' "${dest}" |
  head -n 20 || true

printf 'avx2_ops\t'
grep -E '(vpbroadcast|vpadd|vpcm|vpmov|vpmin|vpmax|vpsad|vpsub|vpcmp|ymm[0-9]+)' "${dest}" |
  head -n 20 || true

printf 'vector_ops\t'
grep -E '(vpadd|vpcm|vpmov|vpmin|vpmax|vpsad|vpsub|vpcmp|movdqa|movdqu|movups|movaps|[xyz]mm[0-9]+)' "${dest}" |
  head -n 20 || true
