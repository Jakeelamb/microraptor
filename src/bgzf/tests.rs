use std::io::{Read, Write};

use super::*;

fn patterned_input(len: usize) -> Vec<u8> {
    (0..len).map(|i| (i % 251) as u8).collect()
}

#[test]
fn bgzf_round_trip_reader_writer() {
    let input = b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nJJJJ\n";
    let mut writer = BgzfWriter::new(Vec::new());
    writer.write_all(input).unwrap();
    let encoded = writer.finish().unwrap();
    assert!(is_bgzf_header(&encoded[..BGZF_HEADER_LEN]));

    let mut reader = BgzfReader::new(&encoded[..]);
    let mut decoded = Vec::new();
    reader.read_to_end(&mut decoded).unwrap();
    assert_eq!(decoded, input);
}

#[test]
fn parallel_round_trip_preserves_order() {
    let mut input = Vec::new();
    for i in 0..10_000 {
        let seq = if i % 2 == 0 { "ACGTACGT" } else { "TGCATGCA" };
        input.extend_from_slice(format!("@r{i}\n{seq}\n+\nIIIIIIII\n").as_bytes());
    }
    let encoded = compress_bgzf_parallel(&input, 4).unwrap();
    let decoded = decompress_bgzf_parallel(&encoded[..], 4).unwrap();
    assert_eq!(decoded, input);
}

#[test]
fn streaming_parallel_reader_round_trip_with_tiny_reads() {
    let mut input = Vec::new();
    for i in 0..1000 {
        input.extend_from_slice(format!("@r{i}\nACGT\n+\nIIII\n").as_bytes());
    }
    let encoded = compress_bgzf_parallel(&input, 4).unwrap();
    let mut reader = BgzfParallelReader::new(std::io::Cursor::new(encoded), 4).unwrap();
    let mut decoded = Vec::new();
    let mut scratch = [0_u8; 37];
    loop {
        let n = reader.read(&mut scratch).unwrap();
        if n == 0 {
            break;
        }
        decoded.extend_from_slice(&scratch[..n]);
    }
    assert_eq!(decoded, input);
}

#[test]
fn reader_rejects_truncated_trailing_header() {
    let input = patterned_input(BGZF_MAX_PAYLOAD + 17);
    let mut encoded = compress_bgzf_parallel(&input, 2).unwrap();
    encoded.truncate(encoded.len() - (BGZF_EOF_BLOCK.len() - 5));

    let mut reader = BgzfReader::new(&encoded[..]);
    let mut decoded = Vec::new();
    let err = reader.read_to_end(&mut decoded).unwrap_err();

    assert_eq!(err.kind(), std::io::ErrorKind::UnexpectedEof);
}

#[test]
fn bounded_send_records_full_queue_metric() {
    let metrics = BgzfPipelineMetrics::default();
    let cancel = AtomicBool::new(false);
    let (tx, rx) = sync_channel(1);
    tx.send(1_u8).unwrap();

    let handle = std::thread::spawn(move || {
        std::thread::sleep(std::time::Duration::from_millis(5));
        assert_eq!(rx.recv().unwrap(), 1);
        assert_eq!(rx.recv().unwrap(), 2);
    });

    assert!(send_bounded(
        &tx,
        2_u8,
        &cancel,
        Some(&metrics),
        BgzfBackpressureChannel::Job,
    ));
    handle.join().unwrap();
    assert!(metrics.snapshot().job_queue_full > 0);
    assert_eq!(metrics.snapshot().result_queue_full, 0);
}

#[test]
fn parallel_reader_accepts_pipeline_metrics() {
    let mut input = Vec::new();
    for i in 0..2000 {
        input.extend_from_slice(format!("@r{i}\nACGT\n+\nIIII\n").as_bytes());
    }
    let encoded = compress_bgzf_parallel(&input, 4).unwrap();
    let metrics = Arc::new(BgzfPipelineMetrics::default());
    let config = BgzfParallelConfig::new(2)
        .with_queue_depths(1, 1)
        .with_metrics(Arc::clone(&metrics));
    let mut reader =
        BgzfParallelReader::with_config(std::io::Cursor::new(encoded), config).unwrap();
    let mut decoded = Vec::new();

    reader.read_to_end(&mut decoded).unwrap();

    assert_eq!(decoded, input);
    let _snapshot = metrics.snapshot();
}

#[test]
fn adaptive_reader_uses_serial_below_parallel_threshold() {
    let input = patterned_input(BGZF_MAX_PAYLOAD + 17);
    let encoded = compress_bgzf_parallel(&input, 2).unwrap();
    let reader = BgzfAutoReader::with_config(
        std::io::Cursor::new(encoded),
        1024,
        BgzfParallelConfig::new(4),
    )
    .unwrap();
    assert!(matches!(reader, BgzfAutoReader::Serial(_)));
}

#[test]
fn adaptive_reader_uses_parallel_above_parallel_threshold() {
    let input = patterned_input(BGZF_MAX_PAYLOAD + 17);
    let encoded = compress_bgzf_parallel(&input, 2).unwrap();
    let reader = BgzfAutoReader::with_config(
        std::io::Cursor::new(encoded),
        DEFAULT_PARALLEL_MIN_COMPRESSED_BYTES,
        BgzfParallelConfig::new(4),
    )
    .unwrap();
    assert!(matches!(reader, BgzfAutoReader::Parallel(_)));
}

#[test]
#[cfg(feature = "libdeflate")]
fn libdeflate_reader_matches_flate2_reader() {
    let mut input = Vec::new();
    for i in 0..2000 {
        input.extend_from_slice(format!("@r{i}\nACGTACGT\n+\nIIIIIIII\n").as_bytes());
    }
    let encoded = compress_bgzf_parallel(&input, 4).unwrap();

    let mut flate2_reader = BgzfReader::new(&encoded[..]);
    let mut flate2_out = Vec::new();
    flate2_reader.read_to_end(&mut flate2_out).unwrap();

    let mut libdeflate_reader =
        BgzfReader::with_inflate_backend(&encoded[..], BgzfInflateBackend::Libdeflate);
    let mut libdeflate_out = Vec::new();
    libdeflate_reader.read_to_end(&mut libdeflate_out).unwrap();

    assert_eq!(flate2_out, input);
    assert_eq!(libdeflate_out, input);
}

#[test]
#[cfg(feature = "libdeflate")]
fn libdeflate_parallel_round_trip_preserves_order() {
    let mut input = Vec::new();
    for i in 0..10_000 {
        input.extend_from_slice(format!("@r{i}\nACGT\n+\nIIII\n").as_bytes());
    }
    let encoded = compress_bgzf_parallel(&input, 4).unwrap();
    let decoded = decompress_bgzf_parallel_with_inflate_backend(
        &encoded[..],
        4,
        BgzfInflateBackend::Libdeflate,
    )
    .unwrap();
    assert_eq!(decoded, input);
}

#[test]
#[cfg(feature = "libdeflate")]
fn libdeflate_writer_round_trip_matches_input() {
    let input = patterned_input(BGZF_MAX_PAYLOAD * 2 + 19);
    let mut writer = BgzfWriter::with_deflate_backend(Vec::new(), BgzfDeflateBackend::Libdeflate);
    writer.write_all(&input).unwrap();
    let encoded = writer.finish().unwrap();

    let mut decoded = Vec::new();
    BgzfReader::new(&encoded[..])
        .read_to_end(&mut decoded)
        .unwrap();
    assert_eq!(decoded, input);
}

#[test]
#[cfg(feature = "libdeflate")]
fn libdeflate_parallel_compress_round_trip_matches_input() {
    let input = patterned_input(BGZF_MAX_PAYLOAD * 3 + 7);
    let encoded =
        compress_bgzf_parallel_with_deflate_backend(&input, 4, BgzfDeflateBackend::Libdeflate)
            .unwrap();
    let decoded = decompress_bgzf_parallel(&encoded[..], 4).unwrap();
    assert_eq!(decoded, input);
}

#[test]
fn builds_bgzf_index_for_multi_block_stream() {
    let input = vec![b'A'; BGZF_MAX_PAYLOAD * 2 + 17];
    let encoded = compress_bgzf_parallel(&input, 3).unwrap();
    let index = build_bgzf_index(&encoded[..]).unwrap();

    assert_eq!(index.len(), 3);
    assert_eq!(index.uncompressed_len(), input.len() as u64);
    assert_eq!(index.compressed_len(), encoded.len() as u64);
    assert_eq!(index.entries()[0].compressed_offset, 0);
    assert_eq!(index.entries()[0].uncompressed_offset, 0);
    assert_eq!(
        index.entries()[1].uncompressed_offset,
        u64::from(index.entries()[0].uncompressed_size)
    );

    let first = index
        .virtual_offset_for_uncompressed_offset(0)
        .unwrap()
        .unwrap();
    assert_eq!(first.compressed_offset(), 0);
    assert_eq!(first.in_block_offset(), 0);

    let second_offset = index.entries()[1].uncompressed_offset + 9;
    let second = index
        .virtual_offset_for_uncompressed_offset(second_offset)
        .unwrap()
        .unwrap();
    assert_eq!(
        second.compressed_offset(),
        index.entries()[1].compressed_offset
    );
    assert_eq!(second.in_block_offset(), 9);
    assert!(
        index
            .virtual_offset_for_uncompressed_offset(input.len() as u64)
            .unwrap()
            .is_none()
    );
}

#[test]
fn builds_empty_index_for_eof_only_stream() {
    let writer = BgzfWriter::new(Vec::new());
    let encoded = writer.finish().unwrap();
    let index = build_bgzf_index(&encoded[..]).unwrap();

    assert!(index.is_empty());
    assert_eq!(index.uncompressed_len(), 0);
    assert_eq!(index.compressed_len(), encoded.len() as u64);
}

#[test]
fn seek_reader_seeks_to_block_start() {
    let input = patterned_input(BGZF_MAX_PAYLOAD * 2 + 123);
    let encoded = compress_bgzf_parallel(&input, 3).unwrap();
    let index = build_bgzf_index(&encoded[..]).unwrap();
    let entry = &index.entries()[1];
    let offset = entry.block_virtual_offset().unwrap();

    let mut reader = BgzfSeekReader::new(std::io::Cursor::new(encoded));
    reader.seek_virtual_offset(offset).unwrap();
    let mut out = vec![0_u8; 257];
    reader.read_exact(&mut out).unwrap();

    let start = entry.uncompressed_offset as usize;
    assert_eq!(out, input[start..start + 257]);
}

#[test]
fn seek_reader_seeks_into_block() {
    let input = patterned_input(BGZF_MAX_PAYLOAD * 2 + 123);
    let encoded = compress_bgzf_parallel(&input, 3).unwrap();
    let index = build_bgzf_index(&encoded[..]).unwrap();
    let entry = &index.entries()[1];
    let in_block = 321;
    let offset = BgzfVirtualOffset::from_parts(entry.compressed_offset, in_block).unwrap();

    let mut reader = BgzfSeekReader::new(std::io::Cursor::new(encoded));
    reader.seek_virtual_offset(offset).unwrap();
    let mut out = vec![0_u8; 4096];
    reader.read_exact(&mut out).unwrap();

    let start = entry.uncompressed_offset as usize + usize::from(in_block);
    assert_eq!(out, input[start..start + 4096]);
}

#[test]
fn seek_reader_rejects_invalid_in_block_offset() {
    let input = patterned_input(BGZF_MAX_PAYLOAD + 1);
    let encoded = compress_bgzf_parallel(&input, 2).unwrap();
    let index = build_bgzf_index(&encoded[..]).unwrap();
    let entry = &index.entries()[0];
    let invalid = u16::try_from(entry.uncompressed_size + 1).unwrap();
    let offset = BgzfVirtualOffset::from_parts(entry.compressed_offset, invalid).unwrap();

    let mut reader = BgzfSeekReader::new(std::io::Cursor::new(encoded));
    let err = reader.seek_virtual_offset(offset).unwrap_err();

    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    assert!(err.to_string().contains("in-block virtual offset"));
}

#[test]
fn seek_reader_reads_to_eof_after_seek() {
    let input = patterned_input(BGZF_MAX_PAYLOAD * 2 + 123);
    let encoded = compress_bgzf_parallel(&input, 3).unwrap();
    let index = build_bgzf_index(&encoded[..]).unwrap();
    let entry = &index.entries()[1];
    let in_block = 17;
    let offset = BgzfVirtualOffset::from_parts(entry.compressed_offset, in_block).unwrap();

    let mut reader = BgzfSeekReader::new(std::io::Cursor::new(encoded));
    reader.seek_virtual_offset(offset).unwrap();
    let mut out = Vec::new();
    reader.read_to_end(&mut out).unwrap();

    let start = entry.uncompressed_offset as usize + usize::from(in_block);
    assert_eq!(out, input[start..]);
}

#[test]
fn virtual_offset_checks_compressed_offset_range() {
    let err = BgzfVirtualOffset::from_parts(1_u64 << 48, 0).unwrap_err();
    assert!(err.to_string().contains("virtual-offset range"));

    let vo = BgzfVirtualOffset::from_parts(123, 45).unwrap();
    assert_eq!(vo.raw(), (123 << 16) | 45);
    assert_eq!(vo.compressed_offset(), 123);
    assert_eq!(vo.in_block_offset(), 45);
}
