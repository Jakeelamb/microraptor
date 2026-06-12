# microraptor benchmark metadata

- generated_at_utc: 2026-06-12T21:49:39Z
- git_commit: fb7166f
- git_dirty: true
- rustc: rustc 1.92.0 (ded5c06cf 2025-12-08)
- cargo_nightly: cargo 1.97.0-nightly (4d1f98451 2026-05-15)
- kernel: Linux 7.0.9-arch2-1 x86_64 GNU/Linux
- cpu: AMD Ryzen AI 9 HX 370 w/ Radeon 890M
- logical_cpus: 24
- memory: 93Gi
- records: 100000
- read_len: 150
- iters: 3
- workers: 24

## Tool versions

### microraptor

```text
microraptor-bench [--input PATH | --paired-inputs R1 R2] [--mode all|parse|pack] [--records N] [--read-len N] [--iters N] [--slab-size BYTES] [--workers N] [--bgzf-parallel-min-bytes N] [--json] [--check-bgzf-pack-regression] [--check-label NAME] [--min-input-bytes N] [--tolerance-pct N] [--skip-timing-checks] [--profile-bgzf-parallel]
```

### seqkit

```text
seqkit v2.13.0
```

### seqtk

```text

Usage:   seqtk <command> <arguments>
Version: 1.5-r133
```

### samtools

```text
samtools 1.23.1
Using htslib 1.23.1
Copyright (C) 2025 Genome Research Ltd.
```

### bgzip

```text
bgzip (htslib) 1.23.1
Copyright (C) 2025 Genome Research Ltd.
```

### fastp

```text
fastp 1.3.3
```
