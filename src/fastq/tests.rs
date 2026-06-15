use super::*;
use crate::FastqPosition;

fn collect_records(input: &[u8], slab_size: usize) -> Result<Vec<(Vec<u8>, Vec<u8>)>> {
    let mut reader = FastqReader::with_config(
        input,
        FastqConfig {
            slab_size,
            validate: true,
            ..FastqConfig::default()
        },
    );
    let mut out = Vec::new();
    while let Some(batch) = reader.next_batch()? {
        for rec in batch.records() {
            out.push((rec.name().to_vec(), rec.seq().to_vec()));
        }
    }
    Ok(out)
}

#[test]
fn reads_single_batch() {
    let out = collect_records(b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nJJJJ\n", 1024).unwrap();
    assert_eq!(
        out,
        vec![
            (b"@r1".to_vec(), b"ACGT".to_vec()),
            (b"@r2".to_vec(), b"TGCA".to_vec())
        ]
    );
}

#[test]
fn handles_crlf() {
    let out = collect_records(b"@r1\r\nACGT\r\n+\r\nIIII\r\n", 1024).unwrap();
    assert_eq!(out, vec![(b"@r1".to_vec(), b"ACGT".to_vec())]);
}

#[test]
fn accepts_missing_final_newline() {
    let out = collect_records(b"@r1\nACGT\n+\nIIII", 1024).unwrap();
    assert_eq!(out, vec![(b"@r1".to_vec(), b"ACGT".to_vec())]);
}

#[test]
fn carries_split_records() {
    let out = collect_records(b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nJJJJ\n", 18).unwrap();
    assert_eq!(out.len(), 2);
    assert_eq!(out[1].1, b"TGCA");
}

#[test]
fn visit_records_matches_batch_records_across_slab_carry() {
    let input = b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nJJJJ\n@r3\nNN\n+\n!!";
    let expected = collect_records(input, 18).unwrap();
    let mut reader = FastqReader::with_config(
        &input[..],
        FastqConfig {
            slab_size: 18,
            validate: true,
            ..FastqConfig::default()
        },
    );
    let mut visited = Vec::new();

    reader
        .visit_records(|record| {
            visited.push((record.name().to_vec(), record.seq().to_vec()));
            Ok(())
        })
        .unwrap();

    assert_eq!(visited, expected);
}

#[test]
fn visit_records_reports_truncated_eof() {
    let mut reader = FastqReader::new(&b"@r1\nACGT\n+"[..]);
    let err = reader.visit_records(|_| Ok(())).unwrap_err();

    assert!(err.to_string().contains("truncated FASTQ record"));
    assert_eq!(error_position(&err), Some(FastqPosition::new(0, 0, 3)));
}

#[test]
fn visit_fastq_bytes_matches_batch_records() {
    let input = b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nJJJJ\n@r3\nNN\n+\n!!";
    let expected = collect_records(input, 18).unwrap();
    let mut visited = Vec::new();

    let records = visit_fastq_bytes(input, FastqConfig::default(), |record| {
        visited.push((record.name().to_vec(), record.seq().to_vec()));
        Ok(())
    })
    .unwrap();

    assert_eq!(records, 3);
    assert_eq!(visited, expected);
}

#[test]
fn visit_fastq_bytes_reports_truncated_eof() {
    let err =
        visit_fastq_bytes(&b"@r1\nACGT\n+"[..], FastqConfig::default(), |_| Ok(())).unwrap_err();

    assert!(err.to_string().contains("truncated FASTQ record"));
    assert_eq!(error_position(&err), Some(FastqPosition::new(0, 0, 3)));
}

#[test]
fn frame_records_carries_partial_line_after_complete_record() {
    let input = b"@r1\nACGT\n+\nIIII\n@partial";
    let mut newlines = Vec::new();
    scan_newlines(input, &mut newlines);
    let mut records = Vec::new();

    let next_start = frame_records(input, &newlines, false, true, 0, 0, &mut records).unwrap();

    assert_eq!(records.len(), 1);
    assert_eq!(next_start, b"@r1\nACGT\n+\nIIII\n".len());
}

#[test]
fn rejects_bad_plus_line() {
    let err = collect_records(b"@r1\nACGT\n-\nIIII\n", 1024).unwrap_err();
    assert!(err.to_string().contains("plus line"));
    assert_eq!(error_position(&err), Some(FastqPosition::new(9, 0, 2)));
}

#[test]
fn rejects_empty_first_token_ids() {
    let err = collect_records(b"@ comment only\nACGT\n+\nIIII\n", 1024).unwrap_err();
    assert!(err.to_string().contains("empty FASTQ id"));
    assert_eq!(error_position(&err), Some(FastqPosition::new(0, 0, 0)));

    let err = visit_fastq_bytes(
        b"@ comment only\nACGT\n+\nIIII\n",
        FastqConfig::default(),
        |_| Ok(()),
    )
    .unwrap_err();
    assert!(err.to_string().contains("empty FASTQ id"));
}

#[test]
fn rejects_quality_length_mismatch() {
    let err = collect_records(b"@r1\nACGT\n+\nIII\n", 1024).unwrap_err();
    assert!(
        err.to_string()
            .contains("quality length 3 != sequence length 4")
    );
    assert_eq!(error_position(&err), Some(FastqPosition::new(11, 0, 3)));
}

#[test]
fn rejects_truncated_eof() {
    let err = collect_records(b"@r1\nACGT\n+", 1024).unwrap_err();
    assert!(err.to_string().contains("truncated FASTQ record"));
    assert_eq!(error_position(&err), Some(FastqPosition::new(0, 0, 3)));
}

#[test]
fn record_id_helpers_are_zero_copy_views() {
    let input = b"@INST:1:FC:2:1101:1000:1000 1:N:0:ATCACG\nACGT\n+\nIIII\n";
    let mut reader = FastqReader::new(&input[..]);
    let batch = reader.next_batch().unwrap().unwrap();
    let rec = batch.records().next().unwrap();

    assert_eq!(rec.name(), b"@INST:1:FC:2:1101:1000:1000 1:N:0:ATCACG");
    assert_eq!(
        rec.name_without_at(),
        b"INST:1:FC:2:1101:1000:1000 1:N:0:ATCACG"
    );
    assert_eq!(rec.id_token(), b"INST:1:FC:2:1101:1000:1000");
    assert_eq!(rec.pair_normalized_id(), b"INST:1:FC:2:1101:1000:1000");
}

#[test]
fn pair_normalized_id_strips_slash_pair_suffixes() {
    let input = b"@frag/1 extra\nACGT\n+\nIIII\n@frag/2 extra\nTGCA\n+\nJJJJ\n";
    let mut reader = FastqReader::new(&input[..]);
    let batch = reader.next_batch().unwrap().unwrap();
    let ids = batch
        .records()
        .map(|r| r.pair_normalized_id().to_vec())
        .collect::<Vec<_>>();
    assert_eq!(ids, vec![b"frag".to_vec(), b"frag".to_vec()]);
}

#[test]
fn strip_pair_suffix_leaves_other_ids_unchanged() {
    assert_eq!(strip_pair_suffix(b"frag/1"), b"frag");
    assert_eq!(strip_pair_suffix(b"frag/2"), b"frag");
    assert_eq!(strip_pair_suffix(b"frag/3"), b"frag/3");
    assert_eq!(strip_pair_suffix(b"frag"), b"frag");
    assert_eq!(strip_pair_suffix(b"1"), b"1");
}

#[test]
fn paired_records_zips_two_batches_by_normalized_id() {
    let r1 = b"@frag1/1\nACGT\n+\nIIII\n@frag2/1\nTGCA\n+\nJJJJ\n";
    let r2 = b"@frag1/2\nACGA\n+\nHHHH\n@frag2/2\nTGCT\n+\nGGGG\n";
    let mut first = FastqReader::new(&r1[..]);
    let mut second = FastqReader::new(&r2[..]);
    let first_batch = first.next_batch().unwrap().unwrap();
    let second_batch = second.next_batch().unwrap().unwrap();

    let pairs = paired_records(&first_batch, &second_batch)
        .unwrap()
        .map(|pair| {
            (
                pair.pair_id().to_vec(),
                pair.first().seq().to_vec(),
                pair.second().seq().to_vec(),
            )
        })
        .collect::<Vec<_>>();

    assert_eq!(
        pairs,
        vec![
            (b"frag1".to_vec(), b"ACGT".to_vec(), b"ACGA".to_vec()),
            (b"frag2".to_vec(), b"TGCA".to_vec(), b"TGCT".to_vec())
        ]
    );
}

#[test]
fn paired_records_rejects_identifier_mismatch() {
    let r1 = b"@frag1/1\nACGT\n+\nIIII\n";
    let r2 = b"@other/2\nACGA\n+\nHHHH\n";
    let mut first = FastqReader::new(&r1[..]);
    let mut second = FastqReader::new(&r2[..]);
    let first_batch = first.next_batch().unwrap().unwrap();
    let second_batch = second.next_batch().unwrap().unwrap();

    let err = paired_records(&first_batch, &second_batch).unwrap_err();
    assert!(err.to_string().contains("identifiers do not match"));
    assert_eq!(error_position(&err), Some(FastqPosition::new(0, 0, 0)));
}

#[test]
fn paired_records_use_first_batch_pair_validation_mode() {
    let r1 = b"@frag1/1\nACGT\n+\nIIII\n";
    let r2 = b"@other/2\nACGA\n+\nHHHH\n";
    let mut first = FastqReader::with_config(
        &r1[..],
        FastqConfig::default().pair_validation(PairValidation::None),
    );
    let mut second = FastqReader::new(&r2[..]);
    let first_batch = first.next_batch().unwrap().unwrap();
    let second_batch = second.next_batch().unwrap().unwrap();

    let pairs = paired_records(&first_batch, &second_batch)
        .unwrap()
        .collect::<Vec<_>>();
    assert_eq!(pairs.len(), 1);
}

#[test]
fn paired_reader_carries_uneven_batch_boundaries() {
    let r1 = [
        make_record("frag1/1", 190),
        make_record("frag2/1", 190),
        make_record("frag3/1", 190),
    ]
    .concat();
    let r2 = [
        make_record("frag1/2", 150),
        make_record("frag2/2", 150),
        make_record("frag3/2", 150),
    ]
    .concat();
    let mut reader = PairedFastqReader::with_configs(
        &r1[..],
        FastqConfig {
            slab_size: 1024,
            ..FastqConfig::default()
        },
        &r2[..],
        FastqConfig {
            slab_size: 1024,
            ..FastqConfig::default()
        },
    );

    let first_ids = {
        let first = reader.next_pair_batch().unwrap().unwrap();
        assert_eq!(first.len(), 2);
        first
            .pairs()
            .map(|pair| pair.pair_id().to_vec())
            .collect::<Vec<_>>()
    };

    let second_ids = {
        let second = reader.next_pair_batch().unwrap().unwrap();
        assert_eq!(second.len(), 1);
        second
            .pairs()
            .map(|pair| pair.pair_id().to_vec())
            .collect::<Vec<_>>()
    };

    assert_eq!(first_ids, vec![b"frag1".to_vec(), b"frag2".to_vec()]);
    assert_eq!(second_ids, vec![b"frag3".to_vec()]);
    assert!(reader.next_pair_batch().unwrap().is_none());
}

#[test]
fn paired_reader_rejects_mismatch_after_retained_boundary() {
    let r1 = [
        make_record("frag1/1", 150),
        make_record("frag2/1", 150),
        make_record("frag3/1", 150),
    ]
    .concat();
    let r2 = [
        make_record("frag1/2", 260),
        make_record("other/2", 260),
        make_record("frag3/2", 260),
    ]
    .concat();
    let mut reader = PairedFastqReader::with_configs(
        &r1[..],
        FastqConfig {
            slab_size: 1024,
            ..FastqConfig::default()
        },
        &r2[..],
        FastqConfig {
            slab_size: 1024,
            ..FastqConfig::default()
        },
    );

    {
        let first = reader.next_pair_batch().unwrap().unwrap();
        assert_eq!(first.len(), 1);
    }

    let err = reader.next_pair_batch().unwrap_err();
    assert!(err.to_string().contains("identifiers do not match"));
    assert_eq!(error_position(&err), Some(FastqPosition::new(533, 1, 0)));
}

#[test]
fn paired_reader_rejects_extra_record_at_eof() {
    let r1 = [make_record("frag1/1", 4), make_record("frag2/1", 4)].concat();
    let r2 = make_record("frag1/2", 4);
    let mut reader = PairedFastqReader::new(&r1[..], &r2[..]);

    {
        let first = reader.next_pair_batch().unwrap().unwrap();
        assert_eq!(first.len(), 1);
    }

    let err = reader.next_pair_batch().unwrap_err();
    assert!(err.to_string().contains("different record counts"));
    assert_eq!(error_position(&err), Some(FastqPosition::new(21, 1, 0)));
}

#[test]
fn paired_reader_reports_truncated_mate_position() {
    let r1 = make_record("frag1/1", 4);
    let r2 = b"@frag1/2\nACGT\n+";
    let mut reader = PairedFastqReader::new(&r1[..], &r2[..]);

    let err = reader.next_pair_batch().unwrap_err();
    assert!(err.to_string().contains("truncated FASTQ record"));
    assert_eq!(error_position(&err), Some(FastqPosition::new(0, 0, 3)));
}

#[test]
fn paired_reader_iterates_successful_pairs() {
    let r1 = [make_record("frag1/1", 4), make_record("frag2/1", 4)].concat();
    let r2 = [make_record("frag1/2", 4), make_record("frag2/2", 4)].concat();
    let mut reader = PairedFastqReader::new(&r1[..], &r2[..]);
    let batch = reader.next_pair_batch().unwrap().unwrap();
    assert!(!batch.is_empty());

    let pairs = batch
        .pairs()
        .map(|pair| {
            (
                pair.pair_id().to_vec(),
                pair.first().seq().to_vec(),
                pair.second().seq().to_vec(),
            )
        })
        .collect::<Vec<_>>();

    assert_eq!(
        pairs,
        vec![
            (b"frag1".to_vec(), b"AAAA".to_vec(), b"AAAA".to_vec()),
            (b"frag2".to_vec(), b"AAAA".to_vec(), b"AAAA".to_vec())
        ]
    );
}

#[test]
fn paired_reader_fast_slash_validation_accepts_ordered_slash_mates() {
    let r1 = [make_record("frag1/1", 4), make_record("frag2/1", 4)].concat();
    let r2 = [make_record("frag1/2", 4), make_record("frag2/2", 4)].concat();
    let mut reader = PairedFastqReader::with_config(
        &r1[..],
        &r2[..],
        FastqConfig::default().pair_validation(PairValidation::FastSlash),
    );
    let batch = reader.next_pair_batch().unwrap().unwrap();
    assert_eq!(batch.len(), 2);
}

#[test]
fn paired_reader_fast_slash_validation_falls_back_for_non_slash_ids() {
    let r1 = b"@frag1\nACGT\n+\nIIII\n";
    let r2 = b"@other\nTGCA\n+\nJJJJ\n";
    let mut reader = PairedFastqReader::with_config(
        &r1[..],
        &r2[..],
        FastqConfig::default().pair_validation(PairValidation::FastSlash),
    );
    let err = reader.next_pair_batch().unwrap_err();
    assert!(err.to_string().contains("identifiers do not match"));
}

#[test]
fn paired_reader_fast_slash_validation_rejects_wrong_mate_suffixes() {
    let r1 = b"@frag1/1\nACGT\n+\nIIII\n";
    let r2 = b"@frag1/1\nTGCA\n+\nJJJJ\n";
    let mut reader = PairedFastqReader::with_config(
        &r1[..],
        &r2[..],
        FastqConfig::default().pair_validation(PairValidation::FastSlash),
    );
    let err = reader.next_pair_batch().unwrap_err();
    assert!(err.to_string().contains("identifiers do not match"));
    assert_eq!(error_position(&err), Some(FastqPosition::new(0, 0, 0)));
}

#[test]
fn paired_reader_can_skip_pair_id_validation() {
    let r1 = b"@frag1/1\nACGT\n+\nIIII\n";
    let r2 = b"@other/2\nTGCA\n+\nJJJJ\n";
    let mut reader = PairedFastqReader::with_config(
        &r1[..],
        &r2[..],
        FastqConfig::default().pair_validation(PairValidation::None),
    );
    let batch = reader.next_pair_batch().unwrap().unwrap();
    assert_eq!(batch.len(), 1);
}

#[test]
fn interleaved_pairs_iterates_valid_pairs() {
    let input = b"@frag1/1\nACGT\n+\nIIII\n@frag1/2\nACGA\n+\nHHHH\n";
    let mut reader = FastqReader::with_config(&input[..], FastqConfig::default().interleaved());
    let batch = reader.next_batch().unwrap().unwrap();

    let pairs = batch
        .interleaved_pairs()
        .unwrap()
        .map(|pair| {
            (
                pair.pair_id().to_vec(),
                pair.first().seq().to_vec(),
                pair.second().seq().to_vec(),
            )
        })
        .collect::<Vec<_>>();

    assert_eq!(
        pairs,
        vec![(b"frag1".to_vec(), b"ACGT".to_vec(), b"ACGA".to_vec())]
    );
}

#[test]
fn interleaved_pairs_rejects_identifier_mismatch() {
    let input = b"@frag1/1\nACGT\n+\nIIII\n@other/2\nACGA\n+\nHHHH\n";
    let mut reader = FastqReader::with_config(&input[..], FastqConfig::default().interleaved());
    let batch = reader.next_batch().unwrap().unwrap();

    let err = batch.interleaved_pairs().unwrap_err();
    assert!(err.to_string().contains("identifiers do not match"));
    assert_eq!(error_position(&err), Some(FastqPosition::new(21, 1, 0)));
}

#[test]
fn interleaved_pairs_use_reader_pair_validation_mode() {
    let input = b"@frag1/1\nACGT\n+\nIIII\n@other/2\nACGA\n+\nHHHH\n";
    let mut reader = FastqReader::with_config(
        &input[..],
        FastqConfig::default()
            .interleaved()
            .pair_validation(PairValidation::None),
    );
    let batch = reader.next_batch().unwrap().unwrap();
    assert_eq!(batch.pair_validation(), PairValidation::None);

    let pairs = batch.interleaved_pairs().unwrap().collect::<Vec<_>>();
    assert_eq!(pairs.len(), 1);
}

#[test]
fn interleaved_fast_slash_validation_rejects_wrong_mate_suffixes() {
    let input = b"@frag1/1\nACGT\n+\nIIII\n@frag1/1\nACGA\n+\nHHHH\n";
    let mut reader = FastqReader::with_config(
        &input[..],
        FastqConfig::default()
            .interleaved()
            .pair_validation(PairValidation::FastSlash),
    );
    let batch = reader.next_batch().unwrap().unwrap();

    let err = batch.interleaved_pairs().unwrap_err();
    assert!(err.to_string().contains("identifiers do not match"));
    assert_eq!(error_position(&err), Some(FastqPosition::new(21, 1, 0)));
}

#[test]
fn interleaved_reader_rejects_odd_eof() {
    let input = b"@frag1/1\nACGT\n+\nIIII\n";
    let mut reader = FastqReader::with_config(&input[..], FastqConfig::default().interleaved());

    let err = reader.next_batch().unwrap_err();
    assert!(err.to_string().contains("unpaired record"));
    assert_eq!(error_position(&err), Some(FastqPosition::new(0, 0, 0)));
}

#[test]
fn interleaved_reader_carries_odd_record_to_next_batch() {
    let input = [
        make_record("frag1/1", 150),
        make_record("frag1/2", 150),
        make_record("frag2/1", 150),
        make_record("frag2/2", 150),
    ]
    .concat();
    let mut reader = FastqReader::with_config(
        &input[..],
        FastqConfig {
            slab_size: 1024,
            ..FastqConfig::default().interleaved()
        },
    );

    let first_ids = {
        let batch = reader.next_batch().unwrap().unwrap();
        assert_eq!(batch.len(), 2);
        batch
            .interleaved_pairs()
            .unwrap()
            .map(|pair| pair.pair_id().to_vec())
            .collect::<Vec<_>>()
    };
    let second_ids = {
        let batch = reader.next_batch().unwrap().unwrap();
        assert_eq!(batch.len(), 2);
        batch
            .interleaved_pairs()
            .unwrap()
            .map(|pair| pair.pair_id().to_vec())
            .collect::<Vec<_>>()
    };

    assert_eq!(first_ids, vec![b"frag1".to_vec()]);
    assert_eq!(second_ids, vec![b"frag2".to_vec()]);
    assert!(reader.next_batch().unwrap().is_none());
}

#[test]
fn reports_absolute_position_after_slab_carry() {
    let input = b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n-\nJJJJ\n";
    let err = collect_records(input, 18).unwrap_err();
    assert_eq!(error_position(&err), Some(FastqPosition::new(25, 1, 2)));
}

fn make_record(name: &str, bases: usize) -> Vec<u8> {
    format!("@{name}\n{}\n+\n{}\n", "A".repeat(bases), "I".repeat(bases)).into_bytes()
}

fn error_position(err: &FastqError) -> Option<FastqPosition> {
    match err {
        FastqError::FormatAt { position, .. } => Some(position.clone()),
        _ => None,
    }
}
