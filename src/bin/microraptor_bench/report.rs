use std::fmt::Write as _;
use std::time::Duration;

#[derive(Debug, Clone, Copy)]
pub(crate) struct ReportConfig<'a> {
    pub(crate) mode: &'a str,
    pub(crate) format: &'a str,
    pub(crate) records: usize,
    pub(crate) read_len: usize,
    pub(crate) iters: usize,
    pub(crate) slab_size: usize,
    pub(crate) workers: usize,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct ReportRow<'a> {
    pub(crate) name: &'a str,
    pub(crate) bytes: usize,
    pub(crate) records: u64,
    pub(crate) bases: u64,
    pub(crate) best: Duration,
    pub(crate) samples: &'a [Duration],
    pub(crate) checksum: u64,
    pub(crate) extras: &'a [(&'static str, u64)],
}

pub(crate) fn print_table(source: &str, input_bytes: usize, rows: &[ReportRow<'_>]) {
    println!("source\t{source}");
    println!("input_bytes\t{input_bytes}");
    println!("name\tinput_mb\tbest_ms\tinput_mib_s\trecords_s\tbases_s\tchecksum");
    for row in rows {
        let secs = row.best.as_secs_f64();
        let mib_s = (row.bytes as f64 / 1_048_576.0) / secs;
        let records_s = row.records as f64 / secs;
        let bases_s = row.bases as f64 / secs;
        println!(
            "{}\t{:.3}\t{:.3}\t{:.3}\t{:.0}\t{:.0}\t{}",
            row.name,
            row.bytes as f64 / 1_048_576.0,
            row.best.as_secs_f64() * 1000.0,
            mib_s,
            records_s,
            bases_s,
            row.checksum
        );
    }
}

pub(crate) fn render_json(
    config: &ReportConfig<'_>,
    source: &str,
    input_bytes: usize,
    rows: &[ReportRow<'_>],
) -> String {
    let mut out = String::new();
    let _ = write!(
        out,
        "{{\"source\":{},\"mode\":{},\"format\":{},\"records\":{},\"read_len\":{},\"iters\":{},\"slab_size\":{},\"workers\":{},\"input_bytes\":{},\"measurements\":[",
        JsonStr(source),
        JsonStr(config.mode),
        JsonStr(config.format),
        config.records,
        config.read_len,
        config.iters,
        config.slab_size,
        config.workers,
        input_bytes
    );
    for (i, row) in rows.iter().enumerate() {
        if i != 0 {
            out.push(',');
        }
        let sample_ns = row
            .samples
            .iter()
            .map(Duration::as_nanos)
            .collect::<Vec<_>>();
        let min_ns = sample_ns
            .iter()
            .copied()
            .min()
            .unwrap_or(row.best.as_nanos());
        let max_ns = sample_ns
            .iter()
            .copied()
            .max()
            .unwrap_or(row.best.as_nanos());
        let median_ns = median_nanos(&sample_ns).unwrap_or(row.best.as_nanos());
        let _ = write!(
            out,
            "{{\"name\":{},\"input_bytes\":{},\"records\":{},\"bases\":{},\"best_ns\":{},\"min_ns\":{},\"median_ns\":{},\"max_ns\":{},\"sample_ns\":[",
            JsonStr(row.name),
            row.bytes,
            row.records,
            row.bases,
            row.best.as_nanos(),
            min_ns,
            median_ns,
            max_ns,
        );
        for (sample_idx, value) in sample_ns.iter().enumerate() {
            if sample_idx != 0 {
                out.push(',');
            }
            let _ = write!(out, "{value}");
        }
        let _ = write!(out, "],\"checksum\":{}", row.checksum);
        for (key, value) in row.extras {
            let _ = write!(out, ",{}:{}", JsonStr(key), value);
        }
        out.push('}');
    }
    out.push_str("]}");
    out
}

fn median_nanos(samples: &[u128]) -> Option<u128> {
    if samples.is_empty() {
        return None;
    }
    let mut sorted = samples.to_vec();
    sorted.sort_unstable();
    let mid = sorted.len() / 2;
    if sorted.len().is_multiple_of(2) {
        Some((sorted[mid - 1] + sorted[mid]) / 2)
    } else {
        Some(sorted[mid])
    }
}

struct JsonStr<'a>(&'a str);

impl std::fmt::Display for JsonStr<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_char('"')?;
        for ch in self.0.chars() {
            match ch {
                '"' => f.write_str("\\\"")?,
                '\\' => f.write_str("\\\\")?,
                '\n' => f.write_str("\\n")?,
                '\r' => f.write_str("\\r")?,
                '\t' => f.write_str("\\t")?,
                ch if ch.is_control() => write!(f, "\\u{:04x}", ch as u32)?,
                ch => f.write_char(ch)?,
            }
        }
        f.write_char('"')
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn json_string_escapes_control_characters() {
        assert_eq!(JsonStr("a\"b\\c\n").to_string(), "\"a\\\"b\\\\c\\n\"");
    }

    #[test]
    fn render_json_includes_mode() {
        let config = ReportConfig {
            mode: "pack",
            format: "fastq",
            records: 0,
            read_len: 0,
            iters: 1,
            slab_size: 1,
            workers: 1,
        };
        let json = render_json(&config, "synthetic", 0, &[]);
        assert!(json.contains("\"mode\":\"pack\""));
    }

    #[test]
    fn render_json_includes_measurement_extras() {
        let config = ReportConfig {
            mode: "all",
            format: "fastq",
            records: 0,
            read_len: 0,
            iters: 1,
            slab_size: 1,
            workers: 1,
        };
        let row = ReportRow {
            name: "bgzf",
            bytes: 1,
            records: 2,
            bases: 3,
            best: Duration::from_nanos(4),
            samples: &[
                Duration::from_nanos(7),
                Duration::from_nanos(4),
                Duration::from_nanos(9),
            ],
            checksum: 5,
            extras: &[("bgzf_job_queue_full", 6)],
        };

        let json = render_json(&config, "synthetic", 1, &[row]);

        assert!(json.contains("\"bgzf_job_queue_full\":6"));
        assert!(json.contains("\"sample_ns\":[7,4,9]"));
        assert!(json.contains("\"median_ns\":7"));
    }
}
