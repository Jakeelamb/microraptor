## Summary

- 

## Evidence

- [ ] `cargo fmt --all -- --check`
- [ ] `RUSTFLAGS="-D warnings" cargo check --lib`
- [ ] `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps`
- [ ] `cargo test --all`
- [ ] `cargo package --allow-dirty`

## Benchmark Or Claim Changes

- [ ] No performance claims changed.
- [ ] Performance claims changed and include exact commands, environment,
      comparator versions, raw JSONL/TSV outputs, rendered summaries, and
      generated figures.

## Public Surface

- [ ] No public API or behavior changed.
- [ ] Public API or behavior changed and README, rustdoc, `CHANGELOG.md`, and
      release docs were updated.
