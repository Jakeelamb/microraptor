#![no_main]

use std::io::Read;

use libfuzzer_sys::fuzz_target;
use microraptor::{
    BgzfParallelReader, BgzfReader, compress_bgzf_parallel, decompress_bgzf_parallel,
};

const MAX_INPUT: usize = 256 * 1024;

fuzz_target!(|data: &[u8]| {
    let input = &data[..data.len().min(MAX_INPUT)];
    let workers = 1 + input.first().copied().unwrap_or(0) as usize % 4;

    let Ok(encoded) = compress_bgzf_parallel(input, workers) else {
        return;
    };

    if let Ok(decoded) = decompress_bgzf_parallel(&encoded[..], workers) {
        assert_eq!(decoded, input);
    }

    let mut serial = BgzfReader::new(&encoded[..]);
    let mut serial_out = Vec::new();
    if serial.read_to_end(&mut serial_out).is_ok() {
        assert_eq!(serial_out, input);
    }

    let Ok(mut parallel) = BgzfParallelReader::new(std::io::Cursor::new(encoded), workers) else {
        return;
    };
    let mut parallel_out = Vec::new();
    if parallel.read_to_end(&mut parallel_out).is_ok() {
        assert_eq!(parallel_out, input);
    }
});
