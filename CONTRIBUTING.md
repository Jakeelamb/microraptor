# Contributing

Microraptor is intended for scientific and systems use. Contributions should
preserve correctness, reproducibility, and honest benchmark boundaries.

## Development Checks

Run focused checks while editing, then run the release gate before publishing:

```bash
cargo fmt --all -- --check
RUSTFLAGS="-D warnings" cargo check --lib
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps
cargo test --all
cargo package --allow-dirty
```

Use nightly only for explicitly nightly surfaces:

```bash
cargo +nightly test --all-features
cargo +nightly clippy --all-targets --all-features -- -D warnings
```

## Benchmark Claims

Do not add or strengthen performance claims without benchmark evidence from the
same commit. A useful benchmark contribution records:

- exact command lines and environment variables;
- hardware, OS, Rust toolchain, and feature flags;
- comparator versions;
- raw JSONL or TSV output;
- generated Markdown/SVG figures from repository scripts;
- checksum, record-count, or base-count parity where relevant.

Machine-local benchmark artifacts are evidence, not universal claims. Keep
parser-only, compression, and workflow-tool comparisons separate.

Use the GitHub benchmark issue template for new public performance wording or
new comparison rows. It asks for the same raw artifacts, command lines, and
claim-boundary checks expected during review.

## Data

Do not commit large biological FASTQ datasets. Use scripts and manifests to
point at local or public datasets, and keep generated large intermediates under
`target/` or another ignored output directory.

## Security And Correctness Reports

Use `SECURITY.md` for vulnerabilities or denial-of-service behavior on
untrusted FASTQ, gzip, or BGZF input. Use the parser bug issue template for
ordinary correctness reports and include a minimal synthetic reproducer instead
of large biological data.
