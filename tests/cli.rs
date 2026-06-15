use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn temp_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "microraptor-cli-{name}-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn write_file(path: &Path, bytes: &[u8]) {
    fs::write(path, bytes).unwrap();
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_microraptor"))
        .args(args)
        .output()
        .unwrap()
}

fn stdout(output: std::process::Output) -> String {
    assert!(
        output.status.success(),
        "status={:?}\nstderr={}",
        output.status.code(),
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

fn stderr_failure(output: std::process::Output) -> String {
    assert!(
        !output.status.success(),
        "expected failure, stdout={}",
        String::from_utf8_lossy(&output.stdout)
    );
    String::from_utf8(output.stderr).unwrap()
}

#[test]
fn stats_reports_fastq_and_fasta_counts() {
    let dir = temp_dir("stats");
    let fastq = dir.join("reads.fastq");
    let fasta = dir.join("refs.fasta");
    write_file(&fastq, b"@r1\nACGT\n+\nIIII\n@r2\nTG\n+\n!!\n");
    write_file(&fasta, b">chr1\nAC\nGT\n>chr2\nTTA\n");

    let fastq_out = stdout(run(&["stats", fastq.to_str().unwrap()]));
    assert!(fastq_out.contains("records\t2\n"));
    assert!(fastq_out.contains("bases\t6\n"));
    assert!(fastq_out.contains("checksum\t"));

    let fasta_out = stdout(run(&[
        "stats",
        "--format",
        "fasta",
        fasta.to_str().unwrap(),
    ]));
    assert!(fasta_out.contains("records\t2\n"));
    assert!(fasta_out.contains("bases\t7\n"));
    assert!(fasta_out.contains("checksum\t"));
}

#[test]
fn checksum_reads_fastq_fasta_and_sam_from_stdin() {
    for (format, input) in [
        ("fastq", b"@r1\nACGT\n+\nIIII\n".as_slice()),
        ("fasta", b">r1\nAC\nGT\n".as_slice()),
        (
            "sam",
            b"@HD\tVN:1.6\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n".as_slice(),
        ),
    ] {
        let mut child = Command::new(env!("CARGO_BIN_EXE_microraptor"))
            .args(["checksum", "--format", format, "-"])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        child.stdin.as_mut().unwrap().write_all(input).unwrap();
        let output = child.wait_with_output().unwrap();
        let out = stdout(output);
        assert!(out.contains("records\t1\n"), "{format}: {out}");
        assert!(out.contains("bases\t4\n"), "{format}: {out}");
    }
}

#[test]
fn fasta_index_and_fetch_round_trip_wrapped_range() {
    let dir = temp_dir("fasta-fetch");
    let fasta = dir.join("refs.fasta");
    let fai = dir.join("refs.fasta.fai");
    write_file(&fasta, b">chr1 desc\nACGT\nTGCA\nAA\n>chr2\nGG\n");

    let index = stdout(run(&["fasta-index", fasta.to_str().unwrap()]));
    assert_eq!(index, "chr1\t10\t11\t4\t5\nchr2\t2\t30\t2\t3\n");
    write_file(&fai, index.as_bytes());

    let fetched = stdout(run(&[
        "fasta-fetch",
        fasta.to_str().unwrap(),
        "--fai",
        fai.to_str().unwrap(),
        "--name",
        "chr1",
        "--start",
        "2",
        "--end",
        "8",
    ]));
    assert_eq!(fetched, "GTTGCA\n");
}

#[test]
fn fasta_partitions_prints_stable_tsv() {
    let dir = temp_dir("fasta-partitions");
    let fasta = dir.join("refs.fasta");
    let fai = dir.join("refs.fasta.fai");
    write_file(&fasta, b">chr1\nACGTTGCAAA\n>chr2\nGG\n");

    let index = stdout(run(&["fasta-index", fasta.to_str().unwrap()]));
    write_file(&fai, index.as_bytes());

    let out = stdout(run(&[
        "fasta-partitions",
        fasta.to_str().unwrap(),
        "--fai",
        fai.to_str().unwrap(),
        "--parts",
        "3",
        "--overlap",
        "1",
    ]));
    assert_eq!(
        out,
        "0\tchr1\t0\t4\t0\t5\t0\n1\tchr1\t4\t8\t3\t9\t1\n2\tchr1\t8\t10\t7\t10\t1\n3\tchr2\t0\t2\t0\t2\t0\n"
    );
}

#[test]
fn fasta_chunks_prints_stable_tsv() {
    let dir = temp_dir("fasta-chunks");
    let fasta = dir.join("refs.fasta");
    let fai = dir.join("refs.fasta.fai");
    write_file(&fasta, b">chr1\nACGT\nTGCA\nAA\n");

    let index = stdout(run(&["fasta-index", fasta.to_str().unwrap()]));
    write_file(&fai, index.as_bytes());

    let out = stdout(run(&[
        "fasta-chunks",
        fasta.to_str().unwrap(),
        "--fai",
        fai.to_str().unwrap(),
        "--name",
        "chr1",
        "--start",
        "2",
        "--end",
        "9",
        "--chunk-bases",
        "3",
    ]));
    assert_eq!(out, "chr1\t2\tGTT\nchr1\t5\tGCA\nchr1\t8\tA\n");
}

#[cfg(feature = "bgzf")]
#[test]
fn verify_bgzf_accepts_valid_bgzf_stream() {
    let dir = temp_dir("verify-bgzf");
    let bgzf = dir.join("reads.fastq.bgz");
    let encoded = microraptor::compress_bgzf_parallel(b"@r1\nACGT\n+\nIIII\n", 1).unwrap();
    write_file(&bgzf, &encoded);

    let out = stdout(run(&["verify-bgzf", bgzf.to_str().unwrap()]));
    assert!(out.contains("status\tok\n"));
    assert!(out.contains("blocks\t"));
}

#[test]
fn cli_rejects_unknown_command() {
    let err = stderr_failure(run(&["not-a-command"]));
    assert!(err.contains("unknown command: not-a-command"));
}

#[test]
fn cli_rejects_missing_required_args_and_format() {
    let err = stderr_failure(run(&["stats"]));
    assert!(err.contains("stats requires a path"));

    let err = stderr_failure(run(&["checksum"]));
    assert!(err.contains("checksum requires --format"));

    let err = stderr_failure(run(&["fasta-index"]));
    assert!(err.contains("fasta-index requires exactly one path argument"));
}

#[test]
fn cli_rejects_unsupported_format() {
    let dir = temp_dir("unsupported-format");
    let input = dir.join("reads.fastq");
    write_file(&input, b"@r1\nACGT\n+\nIIII\n");

    let err = stderr_failure(run(&["stats", "--format", "sam", input.to_str().unwrap()]));
    assert!(err.contains("unsupported stats format: sam"));
}

#[test]
fn cli_rejects_invalid_fetch_args() {
    let dir = temp_dir("invalid-fetch");
    let fasta = dir.join("refs.fasta");
    let fai = dir.join("refs.fasta.fai");
    write_file(&fasta, b">chr1\nACGT\n");
    let index = stdout(run(&["fasta-index", fasta.to_str().unwrap()]));
    write_file(&fai, index.as_bytes());

    let err = stderr_failure(run(&[
        "fasta-fetch",
        fasta.to_str().unwrap(),
        "--fai",
        fai.to_str().unwrap(),
        "--name",
        "chr1",
        "--start",
        "x",
        "--end",
        "2",
    ]));
    assert!(err.contains("--start requires an integer value"));

    let err = stderr_failure(run(&[
        "fasta-fetch",
        fasta.to_str().unwrap(),
        "--fai",
        fai.to_str().unwrap(),
        "--name",
        "chr1",
        "--start",
        "3",
        "--end",
        "2",
    ]));
    assert!(err.contains("FASTA range start must be <= end"));
}

#[cfg(feature = "bgzf")]
#[test]
fn verify_bgzf_rejects_trailing_bytes_after_eof_marker() {
    let dir = temp_dir("verify-bgzf-trailing");
    let bgzf = dir.join("reads.fastq.bgz");
    let mut encoded = microraptor::compress_bgzf_parallel(b"@r1\nACGT\n+\nIIII\n", 1).unwrap();
    encoded.extend_from_slice(b"junk");
    write_file(&bgzf, &encoded);

    let err = stderr_failure(run(&["verify-bgzf", bgzf.to_str().unwrap()]));
    assert!(err.contains("trailing bytes after BGZF EOF marker"));
}

#[cfg(feature = "bgzf")]
#[test]
fn verify_bgzf_rejects_invalid_bgzf() {
    let dir = temp_dir("verify-bgzf-invalid");
    let bgzf = dir.join("reads.fastq.bgz");
    write_file(&bgzf, b"not bgzf");

    let err = stderr_failure(run(&["verify-bgzf", bgzf.to_str().unwrap()]));
    assert!(err.contains("failed to fill whole buffer"));
}
