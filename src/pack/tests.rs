use super::*;

#[test]
fn reports_required_lengths() {
    assert_eq!(packed_base_len(0), 0);
    assert_eq!(packed_base_len(1), 1);
    assert_eq!(packed_base_len(4), 1);
    assert_eq!(packed_base_len(5), 2);

    assert_eq!(bit_mask_len(0), 0);
    assert_eq!(bit_mask_len(1), 1);
    assert_eq!(bit_mask_len(8), 1);
    assert_eq!(bit_mask_len(9), 2);
}

#[test]
fn packs_canonical_bases_four_per_byte() {
    let packed = pack_bases(b"ACGTacgt");
    assert_eq!(packed.bases, vec![0b1110_0100, 0b1110_0100]);
    assert_eq!(packed.n_mask, vec![0]);
    assert_eq!(
        packed.summary,
        BaseSummary {
            len: 8,
            a: 2,
            c: 2,
            g: 2,
            t: 2,
            n: 0,
        }
    );

    let decoded: Vec<_> = (0..packed.len())
        .map(|i| packed_base_at(&packed.bases, &packed.n_mask, i))
        .collect();
    assert_eq!(
        decoded,
        vec![
            Some(PackedBase::A),
            Some(PackedBase::C),
            Some(PackedBase::G),
            Some(PackedBase::T),
            Some(PackedBase::A),
            Some(PackedBase::C),
            Some(PackedBase::G),
            Some(PackedBase::T),
        ]
    );
}

#[test]
fn masks_ambiguous_bases_without_rejecting_record() {
    let packed = pack_bases(b"ANGTXry");
    assert_eq!(packed.bases, vec![0b1110_0000, 0]);
    assert_eq!(packed.n_mask, vec![0b0111_0010]);
    assert_eq!(packed.summary.n, 4);
    assert_eq!(packed.summary.canonical_bases(), 3);

    assert_eq!(
        (0..packed.len())
            .map(|i| packed_base_at(&packed.bases, &packed.n_mask, i))
            .collect::<Vec<_>>(),
        vec![
            Some(PackedBase::A),
            Some(PackedBase::N),
            Some(PackedBase::G),
            Some(PackedBase::T),
            Some(PackedBase::N),
            Some(PackedBase::N),
            Some(PackedBase::N),
        ]
    );
}

#[test]
fn reuses_vec_buffers() {
    let mut bases = Vec::with_capacity(16);
    let mut n_mask = Vec::with_capacity(16);
    bases.extend_from_slice(&[255; 8]);
    n_mask.extend_from_slice(&[255; 8]);

    let summary = pack_bases_into(b"AAAAA", &mut bases, &mut n_mask);
    assert_eq!(summary.len, 5);
    assert_eq!(bases, vec![0, 0]);
    assert_eq!(n_mask, vec![0]);
    assert!(bases.capacity() >= 16);
    assert!(n_mask.capacity() >= 16);
}

#[test]
fn slice_pack_reports_small_buffers() {
    let mut bases = [0; 1];
    let mut n_mask = [0; 1];
    let err = pack_bases_into_slices(b"ACGTA", &mut bases, &mut n_mask).unwrap_err();
    assert_eq!(
        err,
        PackError::OutputTooSmall {
            buffer: PackBuffer::Bases,
            needed: 2,
            provided: 1,
        }
    );
}

#[test]
fn summarizes_phred33_quality() {
    let summary = summarize_qualities(b"!5?I").unwrap();
    assert_eq!(
        summary,
        QualitySummary {
            len: 4,
            min_phred: Some(0),
            max_phred: Some(40),
            sum_phred: 90,
            q20_bases: 3,
            q30_bases: 2,
        }
    );
    assert_eq!(summary.mean_phred(), Some(22.5));
}

#[test]
fn summarizes_long_quality_vector_reduction() {
    let qualities: Vec<u8> = (0..97).map(|i| 33 + (i % 41) as u8).collect();
    let summary = summarize_qualities(&qualities).unwrap();
    let phreds: Vec<u8> = qualities.iter().map(|&byte| byte - 33).collect();

    assert_eq!(summary.len, qualities.len());
    assert_eq!(summary.min_phred, phreds.iter().copied().min());
    assert_eq!(summary.max_phred, phreds.iter().copied().max());
    assert_eq!(
        summary.sum_phred,
        phreds.iter().map(|&phred| u64::from(phred)).sum()
    );
    assert_eq!(
        summary.q20_bases,
        phreds.iter().filter(|&&q| q >= 20).count()
    );
    assert_eq!(
        summary.q30_bases,
        phreds.iter().filter(|&&q| q >= 30).count()
    );
}

#[test]
fn packs_bases_and_summarizes_qualities_together() {
    let mut bases = Vec::new();
    let mut n_mask = Vec::new();
    let summary =
        pack_bases_and_summarize_qualities_into(b"ACGTNN", b"!5?III", &mut bases, &mut n_mask)
            .unwrap();

    assert_eq!(bases, vec![0b1110_0100, 0]);
    assert_eq!(n_mask, vec![0b0011_0000]);
    assert_eq!(summary.bases.len, 6);
    assert_eq!(summary.bases.canonical_bases(), 4);
    assert_eq!(summary.bases.n, 2);
    assert_eq!(
        summary.qualities,
        QualitySummary {
            len: 6,
            min_phred: Some(0),
            max_phred: Some(40),
            sum_phred: 170,
            q20_bases: 5,
            q30_bases: 4,
        }
    );
}

#[test]
fn fused_pack_reused_buffers_do_not_leak_stale_bits() {
    let mut bases = Vec::new();
    let mut n_mask = Vec::new();
    pack_bases_and_summarize_qualities_into(b"NNNNNNNNN", b"IIIIIIIII", &mut bases, &mut n_mask)
        .unwrap();
    assert!(n_mask.iter().any(|&byte| byte != 0));

    let summary =
        pack_bases_and_summarize_qualities_into(b"ACGTA", b"IIIII", &mut bases, &mut n_mask)
            .unwrap();

    assert_eq!(bases, vec![0b1110_0100, 0]);
    assert_eq!(n_mask, vec![0]);
    assert_eq!(summary.bases.n, 0);
    assert_eq!(summary.bases.canonical_bases(), 5);
}

#[test]
fn fused_pack_handles_long_canonical_mixed_case_read() {
    let seq = b"ACGTacgtACGTacgtACGTacgtACGTacgtACGTacgt";
    let qual = vec![b'I'; seq.len()];
    let mut bases = Vec::new();
    let mut n_mask = Vec::new();
    let summary =
        pack_bases_and_summarize_qualities_into(seq, &qual, &mut bases, &mut n_mask).unwrap();

    assert_eq!(summary.bases.len, seq.len());
    assert_eq!(summary.bases.n, 0);
    assert_eq!(summary.bases.a, 10);
    assert_eq!(summary.bases.c, 10);
    assert_eq!(summary.bases.g, 10);
    assert_eq!(summary.bases.t, 10);
    assert!(n_mask.iter().all(|&byte| byte == 0));
    assert_eq!(summary.qualities.sum_phred, 40 * seq.len() as u64);
}

#[test]
fn fused_pack_falls_back_for_different_quality_len() {
    let mut bases = Vec::new();
    let mut n_mask = Vec::new();
    let summary =
        pack_bases_and_summarize_qualities_into(b"ACGT", b"!I", &mut bases, &mut n_mask).unwrap();

    assert_eq!(bases, vec![0b1110_0100]);
    assert_eq!(n_mask, vec![0]);
    assert_eq!(summary.bases.canonical_bases(), 4);
    assert_eq!(
        summary.qualities,
        QualitySummary {
            len: 2,
            min_phred: Some(0),
            max_phred: Some(40),
            sum_phred: 40,
            q20_bases: 1,
            q30_bases: 1,
        }
    );
}

#[test]
fn rejects_non_printable_quality() {
    let err = summarize_qualities(b"I\nI").unwrap_err();
    assert_eq!(
        err,
        PackError::InvalidQuality {
            offset: 1,
            byte: b'\n',
        }
    );
}

#[test]
fn bins_qualities_by_phred_thresholds() {
    let mut bins = Vec::with_capacity(16);
    let summary = bin_qualities_into(b"!+5?I", &[10, 20, 30], &mut bins).unwrap();
    assert_eq!(bins, vec![0, 1, 2, 3, 3]);
    assert_eq!(summary.len, 5);
    assert_eq!(summary.q20_bases, 3);
    assert!(bins.capacity() >= 16);
}

#[test]
fn binning_slice_checks_output_len() {
    let mut bins = [0; 2];
    let err = bin_qualities_into_slice(b"IIII", &[20, 30], &mut bins).unwrap_err();
    assert_eq!(
        err,
        PackError::OutputTooSmall {
            buffer: PackBuffer::QualityBins,
            needed: 4,
            provided: 2,
        }
    );
}

#[test]
fn rejects_unsorted_thresholds() {
    let mut bins = Vec::new();
    let err = bin_qualities_into(b"IIII", &[20, 10], &mut bins).unwrap_err();
    assert_eq!(err, PackError::UnsortedQualityThresholds { index: 1 });
    assert!(bins.is_empty());
}

#[test]
fn trusted_fastq_exposes_packed_buffers_to_sink() {
    let mut seen = Vec::new();
    pack_trusted_fastq(b"@r0\nACGTN\n+\nIIIII\n", |record| {
        seen.push((
            record.name.to_vec(),
            record.bases.to_vec(),
            record.n_mask.to_vec(),
            record.summary,
        ));
        Ok(())
    })
    .unwrap();

    assert_eq!(seen.len(), 1);
    assert_eq!(seen[0].0, b"@r0");
    assert_eq!(seen[0].1, vec![0b1110_0100, 0]);
    assert_eq!(seen[0].2, vec![0b0001_0000]);
    assert_eq!(seen[0].3.bases.n, 1);
}

#[test]
fn trusted_direct_fastq_matches_offset_scan() {
    let input = b"@r0\r\nACGTN\r\n+\r\nIIIII\r\n@r1\nTGCA\n+\n!!!!\n";
    let mut offset = Vec::new();
    let mut direct = Vec::new();

    pack_trusted_fastq(input, |record| {
        offset.push((
            record.name.to_vec(),
            record.bases.to_vec(),
            record.n_mask.to_vec(),
            record.summary,
        ));
        Ok(())
    })
    .unwrap();

    pack_trusted_fastq_direct(input, |record| {
        direct.push((
            record.name.to_vec(),
            record.bases.to_vec(),
            record.n_mask.to_vec(),
            record.summary,
        ));
        Ok(())
    })
    .unwrap();

    assert_eq!(direct, offset);
}

#[test]
fn trusted_direct_stream_handles_slab_carry() {
    let input = b"@r0\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@r1\nNNNN\n+\n!!!!";
    let mut records = 0;
    pack_trusted_fastq_read_direct(
        &input[..],
        FastqConfig {
            slab_size: 16,
            ..FastqConfig::default()
        },
        |record| {
            records += 1;
            assert_eq!(record.summary.bases.len, record.seq.len());
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(records, 2);
}

#[test]
fn trusted_paired_fastq_validates_fast_slash_ids() {
    let r1 = b"@frag/1\nACGT\n+\nIIII\n";
    let r2 = b"@frag/2\nTGCA\n+\nIIII\n";
    let mut pairs = 0;
    pack_trusted_paired_fastq_read(
        &r1[..],
        &r2[..],
        FastqConfig::default(),
        crate::PairValidation::FastSlash,
        |pair| {
            assert_eq!(pair.first.summary.bases.len, 4);
            assert_eq!(pair.second.summary.bases.len, 4);
            pairs += 1;
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(pairs, 1);
}

#[test]
fn trusted_paired_fastq_streams_across_slab_carry() {
    let r1 = b"@frag0/1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@frag1/1\nNNNN\n+\n!!!!\n";
    let r2 = b"@frag0/2\nTGCATGCATGCA\n+\nJJJJJJJJJJJJ\n@frag1/2\nACGT\n+\n####\n";
    let mut pairs = 0;
    pack_trusted_paired_fastq_read(
        &r1[..],
        &r2[..],
        FastqConfig {
            slab_size: 16,
            ..FastqConfig::default()
        },
        crate::PairValidation::FastSlash,
        |pair| {
            assert_eq!(pair.first.summary.bases.len, pair.first.seq.len());
            assert_eq!(pair.second.summary.bases.len, pair.second.seq.len());
            pairs += 1;
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(pairs, 2);
}

#[test]
fn trusted_paired_fastq_rejects_mismatched_ids() {
    let r1 = b"@frag-a/1\nACGT\n+\nIIII\n";
    let r2 = b"@frag-b/2\nTGCA\n+\nIIII\n";
    let err = pack_trusted_paired_fastq_read(
        &r1[..],
        &r2[..],
        FastqConfig::default(),
        crate::PairValidation::FastSlash,
        |_pair| Ok(()),
    )
    .unwrap_err();
    assert!(err.to_string().contains("identifiers do not match"));
}

#[test]
fn trusted_paired_fastq_rejects_different_record_counts() {
    let r1 = b"@frag0/1\nACGT\n+\nIIII\n@frag1/1\nTGCA\n+\nIIII\n";
    let r2 = b"@frag0/2\nTGCA\n+\nIIII\n";
    let err = pack_trusted_paired_fastq_read(
        &r1[..],
        &r2[..],
        FastqConfig {
            slab_size: 16,
            ..FastqConfig::default()
        },
        crate::PairValidation::FastSlash,
        |_pair| Ok(()),
    )
    .unwrap_err();
    assert!(err.to_string().contains("different record counts"));
}

#[test]
fn reports_selected_pack_kernel() {
    #[cfg(feature = "simd")]
    assert!(matches!(
        selected_pack_kernel(),
        PackKernel::PortableSimd | PackKernel::Avx2
    ));
    #[cfg(not(feature = "simd"))]
    assert_eq!(selected_pack_kernel(), PackKernel::Scalar);
}
