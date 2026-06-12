# Security Policy

Microraptor is a FASTQ parsing, packing, and BGZF transport library. Security
reports should focus on behavior that can affect downstream tools processing
untrusted sequencing inputs.

## Supported Versions

The project is pre-1.0. Until the first public release, security fixes target
the current `main` branch.

## Reporting A Vulnerability

Please report suspected vulnerabilities privately through GitHub's private
vulnerability reporting flow when it is enabled for the repository. If that is
not available, open a minimal public issue that states a private report is
needed without including exploit details.

Useful reports include:

- the smallest FASTQ, gzip, or BGZF input that triggers the behavior;
- the exact command or API call used;
- the observed failure mode, such as panic, memory exhaustion, infinite loop,
  data corruption, or incorrect validation;
- feature flags and Rust toolchain;
- whether the input is raw FASTQ, ordinary gzip, or BGZF.

Do not include large biological datasets in reports. Reduce the input to a
minimal synthetic reproducer when possible.

## Scope

In scope:

- parser panics or uncontrolled resource use on malformed inputs;
- incorrect acceptance of malformed four-line FASTQ where validation should
  reject it;
- BGZF decompression, indexing, seeking, or writer bugs that can corrupt data or
  misposition reads;
- pack-path bugs that emit incorrect packed bases, masks, or quality summaries.

Out of scope:

- broad performance claims without a correctness or denial-of-service impact;
- behavior on multiline FASTQ, which is explicitly unsupported;
- reordered paired-end synchronization, which is explicitly unsupported;
- vulnerabilities in external tools used only as benchmark comparators.
