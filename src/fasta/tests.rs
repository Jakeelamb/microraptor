use super::*;

#[test]
fn reads_multiline_fasta_records() {
    let input = b">seq1 description\nACG\nTN\n>seq2\nGG\n";
    let mut reader = FastaReader::new(&input[..]);
    let batch = reader.next_batch().unwrap().unwrap();
    let records: Vec<_> = batch.records().collect();

    assert_eq!(records.len(), 2);
    assert_eq!(records[0].name(), b">seq1 description");
    assert_eq!(records[0].id_token(), b"seq1");
    assert_eq!(records[0].seq(), b"ACGTN");
    assert_eq!(records[1].name_without_gt(), b"seq2");
    assert_eq!(records[1].seq(), b"GG");
    assert!(reader.next_batch().unwrap().is_none());
}

#[test]
fn owned_fasta_batch_can_outlive_reader_batch() {
    let input = b">seq1 description\nACG\nTN\n>seq2\nGG\n";
    let mut reader = FastaReader::new(&input[..]);
    let owned = reader.next_owned_batch().unwrap().unwrap();
    drop(reader);

    let records: Vec<_> = owned
        .records()
        .map(|record| (record.id_token().to_vec(), record.seq().to_vec()))
        .collect();

    assert_eq!(owned.first_record_index(), 0);
    assert_eq!(
        records,
        vec![
            (b"seq1".to_vec(), b"ACGTN".to_vec()),
            (b"seq2".to_vec(), b"GG".to_vec())
        ]
    );
}

#[test]
fn fasta_reader_implements_batch_source_trait() {
    use crate::FastaBatchSource;

    let input = b">seq1\nAC\n";
    let mut reader = FastaReader::new(&input[..]);
    let batch = reader.next_fasta_batch().unwrap().unwrap();
    assert_eq!(batch.records().next().unwrap().seq(), b"AC");
}

#[test]
fn reference_config_matches_reference_opener_tuning() {
    let config = FastaConfig::reference();
    assert_eq!(config.batch_records, 16);
    assert_eq!(config.buffer_size, 256 * 1024);
    assert_eq!(config.expected_seq_len, 1024 * 1024);
}

#[test]
fn carries_header_across_batches() {
    let input = b">seq1\nAC\n>seq2\nGT\n";
    let mut reader = FastaReader::with_config(
        &input[..],
        FastaConfig {
            batch_records: 1,
            ..FastaConfig::default()
        },
    );
    let first = reader.next_batch().unwrap().unwrap();
    let first_record = first.records().next().unwrap();
    assert_eq!(first_record.name_without_gt(), b"seq1");
    assert_eq!(first_record.seq(), b"AC");

    let second = reader.next_batch().unwrap().unwrap();
    let second_record = second.records().next().unwrap();
    assert_eq!(second_record.name_without_gt(), b"seq2");
    assert_eq!(second_record.seq(), b"GT");
    assert_eq!(second.first_record_index(), 1);
}

#[test]
fn rejects_non_header_before_first_record() {
    let mut reader = FastaReader::new(&b"ACGT\n"[..]);
    let err = reader.next_batch().unwrap_err();
    assert!(err.to_string().contains("FASTA record header"));
}

#[test]
fn rejects_empty_header() {
    let mut reader = FastaReader::new(&b">\nACGT\n"[..]);
    let err = reader.next_batch().unwrap_err();
    assert!(err.to_string().contains("empty FASTA id"));
}

#[test]
fn visits_records() {
    let input = b">seq1\nAC\n>seq2\nGT\n";
    let mut reader = FastaReader::new(&input[..]);
    let mut seen = Vec::new();
    reader
        .visit_records(|record| {
            seen.push((record.id_token().to_vec(), record.seq().to_vec()));
            Ok(())
        })
        .unwrap();

    assert_eq!(
        seen,
        vec![
            (b"seq1".to_vec(), b"AC".to_vec()),
            (b"seq2".to_vec(), b"GT".to_vec())
        ]
    );
}

#[test]
fn reports_later_empty_header_index() {
    let input = b">seq1\nAC\n>seq2\nGT\n>\nTT\n";
    let mut reader = FastaReader::with_config(
        &input[..],
        FastaConfig {
            batch_records: 8,
            ..FastaConfig::default()
        },
    );
    let err = reader.next_batch().unwrap_err();
    assert!(err.to_string().contains("record 2"));
}

#[test]
fn visits_resident_fasta_bytes() {
    let input = b"\n>seq1 description\nAC\nGT\n>seq2\nTT\n";
    let mut seen = Vec::new();
    let records = visit_fasta_bytes(input, |record| {
        seen.push((record.id_token().to_vec(), record.seq().to_vec()));
        Ok(())
    })
    .unwrap();

    assert_eq!(records, 2);
    assert_eq!(
        seen,
        vec![
            (b"seq1".to_vec(), b"ACGT".to_vec()),
            (b"seq2".to_vec(), b"TT".to_vec())
        ]
    );
}

#[test]
fn resident_visitor_rejects_non_header_before_first_record() {
    let err = visit_fasta_bytes(b"ACGT\n", |_| Ok(())).unwrap_err();
    assert!(err.to_string().contains("FASTA record header"));
}

#[test]
fn resident_visitor_reports_later_empty_header_index() {
    let input = b">seq1\nAC\n>\nTT\n";
    let err = visit_fasta_bytes(input, |_| Ok(())).unwrap_err();
    assert!(err.to_string().contains("record 1"));
}

#[test]
fn detects_fasta_shape() {
    assert_eq!(detect_fasta_shape(b"").unwrap(), FastaShape::Empty);
    assert_eq!(
        detect_fasta_shape(b">seq1\nAC\n>seq2\nTT\n").unwrap(),
        FastaShape::TwoLine
    );
    assert_eq!(
        detect_fasta_shape(b">seq1\nAC\nGT\n").unwrap(),
        FastaShape::Multiline
    );
    assert_eq!(
        detect_fasta_shape(b">seq1\nAC\n\n>seq2\nTT\n").unwrap(),
        FastaShape::Multiline
    );
    assert!(detect_fasta_shape(b"ACGT\n").is_err());
}

#[test]
fn auto_resident_visitor_matches_robust_multiline_visitor() {
    let input = b">seq1\nAC\nGT\n>seq2\nTT\n";
    let mut auto = Vec::new();
    let mut robust = Vec::new();

    visit_fasta_bytes_auto(input, |record| {
        auto.push((record.id_token().to_vec(), record.seq().to_vec()));
        Ok(())
    })
    .unwrap();
    visit_fasta_bytes(input, |record| {
        robust.push((record.id_token().to_vec(), record.seq().to_vec()));
        Ok(())
    })
    .unwrap();

    assert_eq!(auto, robust);
}

#[test]
fn counts_wrapped_fasta_read_and_bytes() {
    let input = b">seq1\nAC\nGT\n>seq2\nTTA\nA\n";
    let from_read = count_fasta_read(&input[..]).unwrap();
    let from_bytes = count_fasta_bytes(input).unwrap();
    let mut reader = FastaReader::new(&input[..]);
    let from_reader = reader.stats().unwrap();

    assert_eq!(from_read, from_bytes);
    assert_eq!(from_reader, from_read);
    assert_eq!(from_read.records, 2);
    assert_eq!(from_read.bases, 8);
}

#[test]
fn stats_continue_after_partial_batch_without_folding_next_records() {
    let input = b">seq1\nAC\n>seq2\nGT\nTA\n>seq3\nCC\n";
    let mut reader = FastaReader::with_config(
        &input[..],
        FastaConfig {
            batch_records: 1,
            ..FastaConfig::default()
        },
    );
    {
        let first = reader.next_batch().unwrap().unwrap();
        assert_eq!(first.records().next().unwrap().seq(), b"AC");
    }

    let stats = reader.stats().unwrap();
    let mut expected = FastaStats::default();
    expected.observe_sequence(b"GTTA");
    expected.observe_sequence(b"CC");
    assert_eq!(stats, expected);
}

#[test]
fn builds_fasta_index_for_wrapped_reference() {
    let input = b">chr1 description\nACGT\nAC\n>chr2\nTTTT\n";
    let index = build_fasta_index(&input[..]).unwrap();
    assert_eq!(index.len(), 2);

    let chr1 = index.get(b"chr1").unwrap();
    assert_eq!(chr1.len, 6);
    assert_eq!(chr1.offset, 18);
    assert_eq!(chr1.line_bases, 4);
    assert_eq!(chr1.line_width, 5);

    let chr2 = index.get(b"chr2").unwrap();
    assert_eq!(chr2.len, 4);
    assert_eq!(chr2.line_bases, 4);
    assert_eq!(chr2.line_width, 5);

    assert_eq!(
        index.to_fai_string(),
        "chr1\t6\t18\t4\t5\nchr2\t4\t32\t4\t5\n"
    );
}

#[test]
fn fasta_index_rejects_inconsistent_non_final_wrapping() {
    let err = build_fasta_index(&b">chr1\nAC\nACGT\nA\n"[..]).unwrap_err();
    assert!(err.to_string().contains("longer than the first"));
}

#[test]
fn fasta_index_rejects_short_internal_wrapping() {
    let err = build_fasta_index(&b">chr1\nACGT\nAC\nA\n"[..]).unwrap_err();
    assert!(err.to_string().contains("inconsistent wrapping"));
}

#[test]
fn fasta_index_rejects_duplicate_names() {
    let err = build_fasta_index(&b">chr1\nAC\n>chr1 desc\nGT\n"[..]).unwrap_err();
    assert!(err.to_string().contains("duplicate"));
}

#[test]
fn fasta_readers_reject_empty_first_token_ids() {
    let input = b"> description only\nACGT\n";
    let mut reader = FastaReader::new(&input[..]);
    let err = reader.next_batch().unwrap_err();
    assert!(err.to_string().contains("empty FASTA id"));

    let err = visit_fasta_bytes(input, |_| Ok(())).unwrap_err();
    assert!(err.to_string().contains("empty FASTA id"));

    let err = build_fasta_index(&input[..]).unwrap_err();
    assert!(err.to_string().contains("empty FASTA id"));
}

#[test]
#[cfg(feature = "bgzf")]
fn builds_bgzf_aware_fasta_index() {
    let input = b">chr1\nACGT\nAC\n>chr2\nTTTT\n";
    let encoded = crate::compress_bgzf_parallel(input, 2).unwrap();
    let index = build_fasta_index_bgzf(&encoded[..]).unwrap();
    let chr1 = index.get(b"chr1").unwrap();
    assert_eq!(chr1.len, 6);
    let vo = chr1.virtual_offset.unwrap();
    assert_eq!(vo.compressed_offset(), 0);
    assert_eq!(vo.in_block_offset(), 6);
}

#[test]
#[cfg(feature = "bgzf")]
fn bgzf_fasta_index_streams_lines_across_block_boundaries() {
    use std::io::Read;

    let mut input = b">chr1\n".to_vec();
    input.extend(std::iter::repeat_n(b'A', 70_000));
    input.extend_from_slice(b"\n>chr2\nTTTT\n");
    let encoded = crate::compress_bgzf_parallel(&input, 2).unwrap();
    let index = build_fasta_index_bgzf(&encoded[..]).unwrap();

    let chr1 = index.get(b"chr1").unwrap();
    assert_eq!(chr1.len, 70_000);
    assert_eq!(chr1.offset, 6);
    assert_eq!(chr1.line_bases, 70_000);
    assert_eq!(chr1.line_width, 70_001);

    let chr2 = index.get(b"chr2").unwrap();
    assert_eq!(chr2.len, 4);
    let chr2_vo = chr2.virtual_offset.unwrap();
    assert!(chr2_vo.compressed_offset() > 0);

    let mut reader = crate::BgzfSeekReader::new(std::io::Cursor::new(encoded));
    reader.seek_virtual_offset(chr2_vo).unwrap();
    let mut out = [0_u8; 4];
    reader.read_exact(&mut out).unwrap();
    assert_eq!(&out, b"TTTT");
}

#[test]
fn visits_two_line_resident_fasta_bytes() {
    let input = b">seq1 description\nACGT\n>seq2\nTT\n";
    let mut seen = Vec::new();
    let records = visit_two_line_fasta_bytes(input, |record| {
        seen.push((record.id_token().to_vec(), record.seq().to_vec()));
        Ok(())
    })
    .unwrap();

    assert_eq!(records, 2);
    assert_eq!(
        seen,
        vec![
            (b"seq1".to_vec(), b"ACGT".to_vec()),
            (b"seq2".to_vec(), b"TT".to_vec())
        ]
    );
}

#[test]
fn counts_two_line_resident_fasta_bytes() {
    let input = b">seq1\nACGT\n>seq2\nTT\n";
    let stats = count_two_line_fasta_bytes(input).unwrap();
    let mut reference = FastaStats::default();
    visit_two_line_fasta_bytes(input, |record| {
        reference.observe_sequence(record.seq());
        Ok(())
    })
    .unwrap();

    assert_eq!(stats, reference);
    assert_eq!(stats.records, 2);
    assert_eq!(stats.bases, 6);
}

#[test]
fn two_line_visitor_rejects_multiline_fasta() {
    let err = visit_two_line_fasta_bytes(b">seq1\nAC\nGT\n", |_| Ok(())).unwrap_err();
    assert!(err.to_string().contains("followed by a header"));
}

#[test]
fn two_line_resident_counter_rejects_multiline_fasta() {
    let err = count_two_line_fasta_bytes(b">seq1\nAC\nGT\n").unwrap_err();
    assert!(err.to_string().contains("followed by a header"));
}

#[test]
fn visits_two_line_fasta_stream() {
    let input = b">seq1 description\nACGT\n>seq2\nTT\n";
    let mut seen = Vec::new();
    let records = visit_two_line_fasta_read(&input[..], |record| {
        seen.push((record.id_token().to_vec(), record.seq().to_vec()));
        Ok(())
    })
    .unwrap();

    assert_eq!(records, 2);
    assert_eq!(
        seen,
        vec![
            (b"seq1".to_vec(), b"ACGT".to_vec()),
            (b"seq2".to_vec(), b"TT".to_vec())
        ]
    );
}

#[test]
fn counts_two_line_fasta_stream() {
    let input = b">seq1\nACGT\n>seq2\nTT\n";
    let stats = count_two_line_fasta_read(&input[..]).unwrap();
    let mut reference = FastaStats::default();
    visit_two_line_fasta_read(&input[..], |record| {
        reference.observe_sequence(record.seq());
        Ok(())
    })
    .unwrap();

    assert_eq!(stats, reference);
    assert_eq!(stats.records, 2);
    assert_eq!(stats.bases, 6);
}

#[test]
fn two_line_counter_rejects_multiline_fasta() {
    let err = count_two_line_fasta_read(&b">seq1\nAC\nGT\n"[..]).unwrap_err();
    assert!(err.to_string().contains("header must start"));
}

#[test]
fn two_line_stream_rejects_multiline_fasta() {
    let err = visit_two_line_fasta_read(&b">seq1\nAC\nGT\n"[..], |_| Ok(())).unwrap_err();
    assert!(err.to_string().contains("header must start"));
}

#[test]
fn parses_fai_and_fetches_wrapped_range() {
    let input = b">chr1 desc\nACGT\nTGCA\nAA\n>chr2\nGG\n";
    let index = build_fasta_index(&input[..]).unwrap();
    let fai = index.to_fai_string();
    let parsed = FastaIndex::from_fai_str(&fai).unwrap();
    assert_eq!(parsed.to_fai_string(), fai);

    let chr1 = parsed.get(b"chr1").unwrap();
    assert_eq!(chr1.sequence_offset(0).unwrap(), 11);
    assert_eq!(chr1.sequence_offset(4).unwrap(), 16);
    assert_eq!(chr1.sequence_spans(2..8).unwrap(), vec![13..15, 16..20]);

    let mut reader = IndexedFastaReader::new(std::io::Cursor::new(input), parsed);
    assert_eq!(reader.fetch(b"chr1", 2..8).unwrap(), b"GTTGCA");
    assert_eq!(reader.fetch(b"chr1", 10..10).unwrap(), b"");
    assert_eq!(reader.fetch(b"chr2", 0..2).unwrap(), b"GG");
}

#[test]
fn streams_indexed_reference_chunks() {
    let input = b">chr1\nACGT\nTGCA\nAA\n";
    let index = build_fasta_index(&input[..]).unwrap();
    let mut reader = IndexedFastaReader::new(std::io::Cursor::new(input), index);
    let chunks = reader
        .reference_chunks(b"chr1", 2..10, 3)
        .unwrap()
        .collect::<Result<Vec<_>>>()
        .unwrap();

    assert_eq!(
        chunks,
        vec![
            FastaReferenceChunk {
                name: b"chr1".to_vec(),
                global_offset: 2,
                seq: b"GTT".to_vec(),
            },
            FastaReferenceChunk {
                name: b"chr1".to_vec(),
                global_offset: 5,
                seq: b"GCA".to_vec(),
            },
            FastaReferenceChunk {
                name: b"chr1".to_vec(),
                global_offset: 8,
                seq: b"AA".to_vec(),
            },
        ]
    );
}

#[test]
fn streams_indexed_reference_chunks_into_reused_buffer() {
    let input = b">chr1\nACGT\nTGCA\nAA\n";
    let index = build_fasta_index(&input[..]).unwrap();
    let mut reader = IndexedFastaReader::new(std::io::Cursor::new(input), index);
    let mut scratch = Vec::with_capacity(16);
    let scratch_ptr = scratch.as_ptr();
    let mut chunks = Vec::new();

    reader
        .reference_chunks_into(
            b"chr1",
            2..10,
            3,
            &mut scratch,
            &mut |chunk: FastaReferenceChunkRef<'_>| {
                chunks.push((
                    chunk.name.to_vec(),
                    chunk.global_offset,
                    chunk.seq.to_vec(),
                    chunk.seq.as_ptr(),
                ));
                Ok(())
            },
        )
        .unwrap();

    assert_eq!(
        chunks
            .iter()
            .map(|(name, offset, seq, _)| (name.clone(), *offset, seq.clone()))
            .collect::<Vec<_>>(),
        vec![
            (b"chr1".to_vec(), 2, b"GTT".to_vec()),
            (b"chr1".to_vec(), 5, b"GCA".to_vec()),
            (b"chr1".to_vec(), 8, b"AA".to_vec()),
        ]
    );
    assert!(chunks.iter().all(|(_, _, _, ptr)| *ptr == scratch_ptr));
}

#[test]
fn plans_balanced_fasta_partitions_with_overlap() {
    let input = b">chr1\nACGTACGTAC\n>chr2\nTTTTTT\n";
    let index = build_fasta_index(&input[..]).unwrap();
    let partitions = plan_fasta_partitions(&index, FastaPartitionConfig::new(3, 2)).unwrap();

    assert_eq!(
        partitions,
        vec![
            FastaPartition {
                partition_index: 0,
                name: b"chr1".to_vec(),
                core: 0..6,
                fetch: 0..8,
            },
            FastaPartition {
                partition_index: 1,
                name: b"chr1".to_vec(),
                core: 6..10,
                fetch: 4..10,
            },
            FastaPartition {
                partition_index: 2,
                name: b"chr2".to_vec(),
                core: 0..6,
                fetch: 0..6,
            },
        ]
    );
    assert_eq!(partitions[1].core_offset_in_fetch(), 2);
}

#[test]
fn indexed_fetch_handles_line_boundary_and_missing_final_newline() {
    let input = b">chr1\nACGT\nTGCA";
    let index = build_fasta_index(&input[..]).unwrap();
    let mut reader = IndexedFastaReader::new(std::io::Cursor::new(input), index);

    assert_eq!(reader.fetch(b"chr1", 0..4).unwrap(), b"ACGT");
    assert_eq!(reader.fetch(b"chr1", 0..8).unwrap(), b"ACGTTGCA");
    assert_eq!(reader.fetch(b"chr1", 4..8).unwrap(), b"TGCA");
}

#[test]
fn rejects_bad_fai_and_bad_fetch_ranges() {
    assert!(FastaIndex::from_fai_str("chr1\t1\t2\t3\n").is_err());
    assert!(FastaIndex::from_fai_str("chr1\t1\t2\t0\t1\n").is_err());
    assert!(FastaIndex::from_fai_str("chr1\t1\t2\t1\t1\nchr1\t1\t2\t1\t1\n").is_err());

    let input = b">chr1\nACGT\n";
    let index = build_fasta_index(&input[..]).unwrap();
    let mut reader = IndexedFastaReader::new(std::io::Cursor::new(input), index);
    assert!(reader.fetch(b"missing", 0..1).is_err());
    assert!(reader.fetch(b"chr1", Range { start: 3, end: 2 }).is_err());
    assert!(reader.fetch(b"chr1", 0..5).is_err());
}

#[test]
fn rejects_overflowing_fai_offset_math() {
    let index = FastaIndex::from_fai_str(&format!("chr1\t10\t{}\t3\t{}\n", u64::MAX - 1, u64::MAX))
        .unwrap();
    let entry = index.get(b"chr1").unwrap();

    assert!(entry.sequence_offset(3).is_err());
    assert!(entry.sequence_spans(0..4).is_err());
}

#[cfg(feature = "bgzf")]
#[test]
fn fetches_bgzf_fasta_range_using_arbitrary_virtual_offsets() {
    let mut input = b">chr1\n".to_vec();
    let seq = (0..70_010).map(|i| b"ACGT"[i % 4]).collect::<Vec<_>>();
    for chunk in seq.chunks(80) {
        input.extend_from_slice(chunk);
        input.push(b'\n');
    }
    let encoded = crate::compress_bgzf_parallel(&input, 2).unwrap();
    let fasta_index = build_fasta_index_bgzf(&encoded[..]).unwrap();
    let bgzf_index = crate::build_bgzf_index_strict(&encoded[..]).unwrap();
    let mut reader =
        BgzfIndexedFastaReader::new(std::io::Cursor::new(encoded), fasta_index, bgzf_index);

    let fetched = reader.fetch(b"chr1", 69_998..70_006).unwrap();
    assert_eq!(fetched, &seq[69_998..70_006]);
}

#[cfg(feature = "bgzf")]
#[test]
fn bgzf_indexed_fetch_allows_empty_terminal_range_without_indexed_offset() {
    let input = b">chr1\n";
    let encoded = crate::compress_bgzf_parallel(input, 1).unwrap();
    let fasta_index = build_fasta_index_bgzf(&encoded[..]).unwrap();
    let bgzf_index = crate::build_bgzf_index_strict(&encoded[..]).unwrap();

    let mut reader =
        BgzfIndexedFastaReader::new(std::io::Cursor::new(encoded), fasta_index, bgzf_index);
    let mut out = b"stale".to_vec();
    reader.fetch_into(b"chr1", 0..0, &mut out).unwrap();
    assert!(out.is_empty());
}

#[cfg(feature = "bgzf")]
#[test]
fn streams_bgzf_indexed_reference_chunks() {
    let mut input = b">chr1\n".to_vec();
    let seq = (0..70_010).map(|i| b"ACGT"[i % 4]).collect::<Vec<_>>();
    for chunk in seq.chunks(80) {
        input.extend_from_slice(chunk);
        input.push(b'\n');
    }
    let encoded = crate::compress_bgzf_parallel(&input, 2).unwrap();
    let fasta_index = build_fasta_index_bgzf(&encoded[..]).unwrap();
    let bgzf_index = crate::build_bgzf_index_strict(&encoded[..]).unwrap();
    let mut reader =
        BgzfIndexedFastaReader::new(std::io::Cursor::new(encoded), fasta_index, bgzf_index);
    let chunks = reader
        .reference_chunks(b"chr1", 69_998..70_006, 5)
        .unwrap()
        .collect::<Result<Vec<_>>>()
        .unwrap();

    assert_eq!(chunks.len(), 2);
    assert_eq!(chunks[0].global_offset, 69_998);
    assert_eq!(chunks[0].seq, &seq[69_998..70_003]);
    assert_eq!(chunks[1].global_offset, 70_003);
    assert_eq!(chunks[1].seq, &seq[70_003..70_006]);
}

#[cfg(feature = "bgzf")]
#[test]
fn streams_bgzf_indexed_reference_chunks_into_reused_buffer() {
    let mut input = b">chr1\n".to_vec();
    let seq = (0..70_010).map(|i| b"ACGT"[i % 4]).collect::<Vec<_>>();
    for chunk in seq.chunks(80) {
        input.extend_from_slice(chunk);
        input.push(b'\n');
    }
    let encoded = crate::compress_bgzf_parallel(&input, 2).unwrap();
    let fasta_index = build_fasta_index_bgzf(&encoded[..]).unwrap();
    let bgzf_index = crate::build_bgzf_index_strict(&encoded[..]).unwrap();
    let mut reader =
        BgzfIndexedFastaReader::new(std::io::Cursor::new(encoded), fasta_index, bgzf_index);
    let mut scratch = Vec::with_capacity(8);
    let scratch_ptr = scratch.as_ptr();
    let mut chunks = Vec::new();

    reader
        .reference_chunks_into(
            b"chr1",
            69_998..70_006,
            5,
            &mut scratch,
            &mut |chunk: FastaReferenceChunkRef<'_>| {
                chunks.push((chunk.global_offset, chunk.seq.to_vec(), chunk.seq.as_ptr()));
                Ok(())
            },
        )
        .unwrap();

    assert_eq!(chunks.len(), 2);
    assert_eq!(chunks[0].0, 69_998);
    assert_eq!(chunks[0].1, &seq[69_998..70_003]);
    assert_eq!(chunks[1].0, 70_003);
    assert_eq!(chunks[1].1, &seq[70_003..70_006]);
    assert!(chunks.iter().all(|(_, _, ptr)| *ptr == scratch_ptr));
}
