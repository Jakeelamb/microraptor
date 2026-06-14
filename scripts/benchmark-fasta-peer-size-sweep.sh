#!/usr/bin/env bash
set -euo pipefail

out_dir="${MICRORAPTOR_FASTA_SIZE_SWEEP_OUT_DIR:-target/bench-results/fasta-peer-size-sweep}"
sizes="${MICRORAPTOR_FASTA_SIZE_SWEEP_RECORDS:-10000 50000 100000 500000 1000000}"
read_len="${MICRORAPTOR_FASTA_SIZE_SWEEP_READ_LEN:-150}"
iters="${MICRORAPTOR_FASTA_SIZE_SWEEP_ITERS:-3}"
consumer="${MICRORAPTOR_FASTA_SIZE_SWEEP_CONSUMER:-light}"
compressions="${MICRORAPTOR_FASTA_SIZE_SWEEP_COMPRESSIONS:-raw gzip}"
base_project_dir="${MICRORAPTOR_FASTA_SIZE_SWEEP_PROJECT_DIR:-target/fasta-peer-size-sweep}"
summary="${out_dir}/summary.md"
tsv="${out_dir}/fasta-peer-size-sweep.tsv"
metadata="${out_dir}/metadata.md"
fig_dir="${out_dir}/figures"
svg="${fig_dir}/fasta-peer-size-sweep-time.svg"
bench_threads="${MICRORAPTOR_BENCH_THREADS:-8}"
require_microraptor_wins="${MICRORAPTOR_FASTA_SIZE_SWEEP_REQUIRE_MICRORAPTOR_WINS:-1}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/benchmark-common.sh"
microraptor_set_thread_cap "${bench_threads}"

mkdir -p "${out_dir}" "${fig_dir}"

tmp_tsv="$(mktemp "${out_dir}/fasta-peer-size-sweep.raw.XXXXXX")"
trap 'rm -f "${tmp_tsv}"' EXIT

printf 'compression\trecords\tread_len\ttool\tbases\tinput_bytes\tbest_ms\trecords_s\tbases_s\tchecksum\titers\tsource\n' > "${tmp_tsv}"

for compression in ${compressions}; do
  for records in ${sizes}; do
    run_dir="${out_dir}/runs/${compression}-${records}"
    project_dir="${base_project_dir}"
    mkdir -p "${run_dir}"
    printf 'running FASTA peer size sweep: compression=%s records=%s read_len=%s iters=%s threads=%s\n' \
      "${compression}" "${records}" "${read_len}" "${iters}" "${bench_threads}"
    MICRORAPTOR_FASTA_PEER_OUT_DIR="${run_dir}" \
    MICRORAPTOR_FASTA_PEER_PROJECT_DIR="${project_dir}" \
    MICRORAPTOR_FASTA_PEER_RECORDS="${records}" \
    MICRORAPTOR_FASTA_PEER_READ_LEN="${read_len}" \
    MICRORAPTOR_FASTA_PEER_ITERS="${iters}" \
    MICRORAPTOR_FASTA_PEER_CONSUMER="${consumer}" \
    MICRORAPTOR_FASTA_PEER_COMPRESSION="${compression}" \
    MICRORAPTOR_BENCH_THREADS="${bench_threads}" \
      scripts/benchmark-fasta-peers.sh >/dev/null

    awk -F '\t' -v compression="${compression}" -v records="${records}" -v read_len="${read_len}" '
      NR > 1 {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
          compression, records, read_len, $1, $3, $8, $4, $5, $6, $7, $9, $10
      }
    ' "${run_dir}/fasta-library-peers.tsv" >> "${tmp_tsv}"
  done
done

mv "${tmp_tsv}" "${tsv}"
trap - EXIT

if [[ "${require_microraptor_wins}" != "0" ]]; then
  awk -F '\t' '
    NR == 1 { next }
    {
      key = $1 "\t" $2
      if (!(key in best) || $7 + 0 < best[key]) {
        best[key] = $7 + 0
        tool[key] = $4
        row[key] = $0
      }
    }
    END {
      failed = 0
      for (key in tool) {
        if (tool[key] !~ /^microraptor/) {
          print "non-microraptor FASTA sweep winner: " row[key] > "/dev/stderr"
          failed = 1
        }
      }
      exit failed
    }
  ' "${tsv}"
fi

awk -F '\t' '
  NR > 1 {
    key = $1 "/" $4
    if (!(key in seen_tool)) {
      tool[++tool_count] = key
      seen_tool[key] = tool_count
    }
    point_count[key]++
    x[key, point_count[key]] = $6 + 0
    y[key, point_count[key]] = $7 + 0
    records[key, point_count[key]] = $2
    if ($6 + 0 > max_x) max_x = $6 + 0
    if ($7 + 0 > max_y) max_y = $7 + 0
  }
  END {
    width = 980
    height = 620
    left = 88
    right = 34
    top = 70
    bottom = 86
    plot_w = width - left - right
    plot_h = height - top - bottom
    split("#1f5a85 #00866a #8a4f00 #7a3f87 #4d5f00 #a23b3b", palette, " ")

    print "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" width "\" height=\"" height "\" viewBox=\"0 0 " width " " height "\">"
    print "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>"
    print "<text x=\"28\" y=\"32\" font-family=\"Arial, sans-serif\" font-size=\"20\" font-weight=\"700\">FASTA parser framework size sweep</text>"
    print "<text x=\"28\" y=\"54\" font-family=\"Arial, sans-serif\" font-size=\"12\" fill=\"#555\">Best wall time over synthetic FASTA inputs; lower is better</text>"
    print "<line x1=\"" left "\" y1=\"" (top + plot_h) "\" x2=\"" (left + plot_w) "\" y2=\"" (top + plot_h) "\" stroke=\"#222\"/>"
    print "<line x1=\"" left "\" y1=\"" top "\" x2=\"" left "\" y2=\"" (top + plot_h) "\" stroke=\"#222\"/>"

    for (i = 0; i <= 4; i++) {
      gy = top + plot_h - (plot_h * i / 4)
      value = max_y * i / 4
      printf "<line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" stroke=\"#e6e6e6\"/>\n", left, gy, left + plot_w, gy
      printf "<text x=\"%d\" y=\"%.1f\" font-family=\"Arial, sans-serif\" font-size=\"11\" text-anchor=\"end\" fill=\"#555\">%.0f ms</text>\n", left - 8, gy + 4, value
    }

    for (i = 0; i <= 4; i++) {
      gx = left + plot_w * i / 4
      value = max_x * i / 4 / 1000000.0
      printf "<line x1=\"%.1f\" y1=\"%d\" x2=\"%.1f\" y2=\"%d\" stroke=\"#e6e6e6\"/>\n", gx, top, gx, top + plot_h
      printf "<text x=\"%.1f\" y=\"%d\" font-family=\"Arial, sans-serif\" font-size=\"11\" text-anchor=\"middle\" fill=\"#555\">%.0f MB</text>\n", gx, top + plot_h + 22, value
    }

    for (t = 1; t <= tool_count; t++) {
      name = tool[t]
      color = palette[((t - 1) % 6) + 1]
      path = ""
      for (i = 1; i <= point_count[name]; i++) {
        px = left + (max_x > 0 ? x[name, i] / max_x * plot_w : 0)
        py = top + plot_h - (max_y > 0 ? y[name, i] / max_y * plot_h : 0)
        path = path sprintf("%s%.1f %.1f", i == 1 ? "M " : " L ", px, py)
      }
      printf "<path d=\"%s\" fill=\"none\" stroke=\"%s\" stroke-width=\"2.5\"/>\n", path, color
      for (i = 1; i <= point_count[name]; i++) {
        px = left + (max_x > 0 ? x[name, i] / max_x * plot_w : 0)
        py = top + plot_h - (max_y > 0 ? y[name, i] / max_y * plot_h : 0)
        printf "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"4\" fill=\"%s\"/>\n", px, py, color
      }
      ly = top + 18 + (t - 1) * 22
      printf "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"%s\" stroke-width=\"3\"/>\n", left + plot_w - 220, ly - 4, left + plot_w - 196, ly - 4, color
      printf "<text x=\"%d\" y=\"%d\" font-family=\"Arial, sans-serif\" font-size=\"12\" fill=\"#222\">%s</text>\n", left + plot_w - 188, ly, name
    }

    print "<text x=\"" (left + plot_w / 2) "\" y=\"" (height - 25) "\" font-family=\"Arial, sans-serif\" font-size=\"13\" text-anchor=\"middle\" fill=\"#222\">Input size, resident FASTA bytes</text>"
    print "<text x=\"22\" y=\"" (top + plot_h / 2) "\" font-family=\"Arial, sans-serif\" font-size=\"13\" text-anchor=\"middle\" fill=\"#222\" transform=\"rotate(-90 22 " (top + plot_h / 2) ")\">Best wall time, ms</text>"
    print "</svg>"
  }
' "${tsv}" > "${svg}"

{
  printf '# FASTA Rust Peer Size Sweep\n\n'
  printf 'Generated by `%s`.\n\n' "$0"
  printf 'This sweep compares Rust FASTA parser frameworks on deterministic synthetic two-line FASTA inputs of increasing size. Raw rows measure one resident byte buffer per size. Gzip rows measure one gzip-compressed resident byte buffer per size; rows explicitly name whether they use third-party flate2 streaming, third-party flate2 resident decode, third-party libdeflate resident decode through the libdeflater Rust wrapper, strict two-line record views, or the light-accounting two-line counter path. This is parser-framework evidence, not command-line workflow evidence for indexing, filtering, BGZF, or a new microraptor DEFLATE implementation.\n\n'
  printf -- '- sizes_records: `%s`\n' "${sizes}"
  printf -- '- read_len: `%s`\n' "${read_len}"
  printf -- '- iters: `%s`\n' "${iters}"
  printf -- '- consumer: `%s`\n' "${consumer}"
  printf -- '- compressions: `%s`\n' "${compressions}"
  printf -- '- thread_cap: `%s`\n\n' "${bench_threads}"
  printf -- '- require_microraptor_wins: `%s`\n\n' "${require_microraptor_wins}"
  printf 'Figure: [`figures/fasta-peer-size-sweep-time.svg`](figures/fasta-peer-size-sweep-time.svg)\n\n'
  printf '| compression | records | input MB | tool | best ms | records/s | bases/s |\n'
  printf '| --- | ---: | ---: | --- | ---: | ---: | ---: |\n'
  awk -F '\t' 'NR > 1 {
    printf "| `%s` | %s | %.2f | `%s` | %.3f | %s | %s |\n", $1, $2, $6 / 1000000.0, $4, $7, $8, $9
  }' "${tsv}"
  printf '\nRaw data: [`fasta-peer-size-sweep.tsv`](fasta-peer-size-sweep.tsv)\n'
  printf 'Metadata: [`metadata.md`](metadata.md)\n'
} > "${summary}"

{
  printf '# FASTA Rust Peer Size Sweep Metadata\n\n'
  printf -- '- generated_at_utc: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf -- '- git_commit: %s\n' "$(git rev-parse --short HEAD 2>/dev/null || printf unknown)"
  printf -- '- git_dirty: %s\n' "$(if [[ -n "$(git status --short 2>/dev/null)" ]]; then printf true; else printf false; fi)"
  printf -- '- rustc: %s\n' "$(rustc --version)"
  printf -- '- cargo: %s\n' "$(cargo --version)"
  printf -- '- kernel: %s\n' "$(uname -srmo)"
  printf -- '- cpu: %s\n' "$(lscpu | sed -n 's/^Model name:[[:space:]]*//p' | head -n 1)"
  printf -- '- logical_cpus: %s\n' "$(nproc)"
  printf -- '- memory: %s\n' "$(free -h | awk '/^Mem:/ { print $2 }')"
  printf -- '- sizes_records: %s\n' "${sizes}"
  printf -- '- read_len: %s\n' "${read_len}"
  printf -- '- iters: %s\n' "${iters}"
  printf -- '- consumer: %s\n' "${consumer}"
  printf -- '- compressions: %s\n' "${compressions}"
  printf -- '- thread_cap: %s\n' "${bench_threads}"
  printf -- '- require_microraptor_wins: %s\n' "${require_microraptor_wins}"
} > "${metadata}"

printf 'wrote %s\n' "${tsv}"
printf 'wrote %s\n' "${summary}"
printf 'wrote %s\n' "${svg}"
printf 'wrote %s\n' "${metadata}"
