use std::io::{Read, Write};

use crc32fast::Hasher;
use flate2::Compression;
use flate2::read::DeflateDecoder;
use flate2::write::DeflateEncoder;

use crate::error::{FastqError, Result};

const BGZF_HEADER_LEN: usize = 18;
const GZIP_TRAILER_LEN: usize = 8;
const BGZF_MAX_BLOCK_SIZE: usize = 64 * 1024;
const BGZF_MAX_PAYLOAD: usize = 60 * 1024;

pub const BGZF_EOF_BLOCK: &[u8] = &[
    31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 27, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

#[derive(Debug)]
struct CompressedBlock {
    bytes: Vec<u8>,
}

impl CompressedBlock {
    fn is_eof(&self) -> bool {
        self.bytes == BGZF_EOF_BLOCK
    }
}

pub fn is_bgzf_header(prefix: &[u8]) -> bool {
    prefix.len() >= BGZF_HEADER_LEN
        && prefix[0] == 31
        && prefix[1] == 139
        && prefix[2] == 8
        && prefix[3] & 4 != 0
        && u16::from_le_bytes([prefix[10], prefix[11]]) >= 6
        && prefix[12] == b'B'
        && prefix[13] == b'C'
        && prefix[14] == 2
        && prefix[15] == 0
}

pub struct BgzfReader<R> {
    inner: R,
    decoded: Vec<u8>,
    pos: usize,
    eof: bool,
}

impl<R: Read> BgzfReader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            decoded: Vec::new(),
            pos: 0,
            eof: false,
        }
    }

    fn refill(&mut self) -> std::io::Result<()> {
        self.decoded.clear();
        self.pos = 0;
        loop {
            let Some(block) = read_block(&mut self.inner)? else {
                self.eof = true;
                return Ok(());
            };
            if block.is_eof() {
                self.eof = true;
                return Ok(());
            }
            self.decoded = decode_block(&block)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            if !self.decoded.is_empty() {
                return Ok(());
            }
        }
    }
}

impl<R: Read> Read for BgzfReader<R> {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        if self.pos >= self.decoded.len() {
            if self.eof {
                return Ok(0);
            }
            self.refill()?;
            if self.pos >= self.decoded.len() && self.eof {
                return Ok(0);
            }
        }
        let n = out.len().min(self.decoded.len() - self.pos);
        out[..n].copy_from_slice(&self.decoded[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

pub struct BgzfWriter<W> {
    inner: Option<W>,
    pending: Vec<u8>,
    level: Compression,
}

impl<W: Write> BgzfWriter<W> {
    pub fn new(inner: W) -> Self {
        Self::with_compression(inner, Compression::fast())
    }

    pub fn with_compression(inner: W, level: Compression) -> Self {
        Self {
            inner: Some(inner),
            pending: Vec::with_capacity(BGZF_MAX_PAYLOAD),
            level,
        }
    }

    pub fn finish(mut self) -> Result<W> {
        self.flush_pending()?;
        let mut inner = self
            .inner
            .take()
            .ok_or_else(|| FastqError::Bgzf("writer already finished".into()))?;
        inner.write_all(BGZF_EOF_BLOCK)?;
        inner.flush()?;
        Ok(inner)
    }

    fn inner_mut(&mut self) -> std::io::Result<&mut W> {
        self.inner
            .as_mut()
            .ok_or_else(|| std::io::Error::other("BGZF writer already finished"))
    }

    fn flush_pending(&mut self) -> Result<()> {
        if self.pending.is_empty() {
            return Ok(());
        }
        let block = encode_block(&self.pending, self.level)?;
        self.inner_mut()?.write_all(&block)?;
        self.pending.clear();
        Ok(())
    }
}

impl<W: Write> Write for BgzfWriter<W> {
    fn write(&mut self, mut buf: &[u8]) -> std::io::Result<usize> {
        let original = buf.len();
        while !buf.is_empty() {
            let space = BGZF_MAX_PAYLOAD - self.pending.len();
            let take = space.min(buf.len());
            self.pending.extend_from_slice(&buf[..take]);
            buf = &buf[take..];
            if self.pending.len() == BGZF_MAX_PAYLOAD {
                self.flush_pending().map_err(std::io::Error::other)?;
            }
        }
        Ok(original)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.flush_pending().map_err(std::io::Error::other)?;
        self.inner_mut()?.flush()
    }
}

pub fn decompress_bgzf_parallel<R: Read>(mut reader: R, workers: usize) -> Result<Vec<u8>> {
    let mut blocks = Vec::new();
    while let Some(block) = read_block(&mut reader)? {
        if block.is_eof() {
            break;
        }
        blocks.push(block);
    }
    let decoded = parallel_map(blocks, workers, |block| decode_block(&block))?;
    let total = decoded.iter().map(Vec::len).sum();
    let mut out = Vec::with_capacity(total);
    for chunk in decoded {
        out.extend_from_slice(&chunk);
    }
    Ok(out)
}

pub fn compress_bgzf_parallel(input: &[u8], workers: usize) -> Result<Vec<u8>> {
    let chunks: Vec<&[u8]> = input.chunks(BGZF_MAX_PAYLOAD).collect();
    let encoded = parallel_map(chunks, workers, |chunk| {
        encode_block(chunk, Compression::fast())
    })?;
    let total = encoded.iter().map(Vec::len).sum::<usize>() + BGZF_EOF_BLOCK.len();
    let mut out = Vec::with_capacity(total);
    for block in encoded {
        out.extend_from_slice(&block);
    }
    out.extend_from_slice(BGZF_EOF_BLOCK);
    Ok(out)
}

fn parallel_map<T, F>(items: Vec<T>, workers: usize, f: F) -> Result<Vec<Vec<u8>>>
where
    T: Send,
    F: Fn(T) -> Result<Vec<u8>> + Sync,
{
    let n = items.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let worker_count = workers.max(1).min(n);
    let mut buckets = (0..worker_count)
        .map(|_| Vec::new())
        .collect::<Vec<Vec<(usize, T)>>>();
    for (idx, item) in items.into_iter().enumerate() {
        buckets[idx % worker_count].push((idx, item));
    }

    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for bucket in buckets {
            let f_ref = &f;
            handles.push(scope.spawn(move || {
                let mut out = Vec::with_capacity(bucket.len());
                for (idx, item) in bucket {
                    out.push((idx, f_ref(item)?));
                }
                Result::<Vec<(usize, Vec<u8>)>>::Ok(out)
            }));
        }

        let mut ordered = vec![Vec::new(); n];
        for handle in handles {
            let rows = handle
                .join()
                .map_err(|_| FastqError::Bgzf("parallel worker panicked".into()))??;
            for (idx, bytes) in rows {
                ordered[idx] = bytes;
            }
        }
        Ok(ordered)
    })
}

fn read_block<R: Read>(reader: &mut R) -> std::io::Result<Option<CompressedBlock>> {
    let mut header = [0_u8; BGZF_HEADER_LEN];
    match reader.read_exact(&mut header) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e),
    }
    if !is_bgzf_header(&header) {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "invalid BGZF header",
        ));
    }
    let bsize = u16::from_le_bytes([header[16], header[17]]) as usize + 1;
    if !(BGZF_HEADER_LEN + GZIP_TRAILER_LEN..=BGZF_MAX_BLOCK_SIZE).contains(&bsize) {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "invalid BGZF block size",
        ));
    }
    let mut bytes = Vec::with_capacity(bsize);
    bytes.extend_from_slice(&header);
    bytes.resize(bsize, 0);
    reader.read_exact(&mut bytes[BGZF_HEADER_LEN..])?;
    Ok(Some(CompressedBlock { bytes }))
}

fn decode_block(block: &CompressedBlock) -> Result<Vec<u8>> {
    let bytes = &block.bytes;
    if bytes.len() < BGZF_HEADER_LEN + GZIP_TRAILER_LEN {
        return Err(FastqError::Bgzf("short block".into()));
    }
    let compressed_end = bytes.len() - GZIP_TRAILER_LEN;
    let mut decoder = DeflateDecoder::new(&bytes[BGZF_HEADER_LEN..compressed_end]);
    let mut out = Vec::new();
    decoder.read_to_end(&mut out)?;

    let expected_crc = u32::from_le_bytes([
        bytes[compressed_end],
        bytes[compressed_end + 1],
        bytes[compressed_end + 2],
        bytes[compressed_end + 3],
    ]);
    let expected_len = u32::from_le_bytes([
        bytes[compressed_end + 4],
        bytes[compressed_end + 5],
        bytes[compressed_end + 6],
        bytes[compressed_end + 7],
    ]) as usize;
    if out.len() != expected_len {
        return Err(FastqError::Bgzf("uncompressed size mismatch".into()));
    }
    let mut hasher = Hasher::new();
    hasher.update(&out);
    if hasher.finalize() != expected_crc {
        return Err(FastqError::Bgzf("CRC32 mismatch".into()));
    }
    Ok(out)
}

fn encode_block(input: &[u8], level: Compression) -> Result<Vec<u8>> {
    if input.len() > BGZF_MAX_PAYLOAD {
        return Err(FastqError::Bgzf("payload exceeds BGZF block policy".into()));
    }

    let mut encoder = DeflateEncoder::new(Vec::new(), level);
    encoder.write_all(input)?;
    let compressed = encoder.finish()?;
    let total_size = BGZF_HEADER_LEN + compressed.len() + GZIP_TRAILER_LEN;
    if total_size > BGZF_MAX_BLOCK_SIZE {
        return Err(FastqError::Bgzf(
            "compressed BGZF block exceeds 64 KiB".into(),
        ));
    }

    let mut out = Vec::with_capacity(total_size);
    out.extend_from_slice(&[31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0]);
    out.extend_from_slice(&[b'B', b'C', 2, 0]);
    out.extend_from_slice(
        &u16::try_from(total_size - 1)
            .unwrap_or(u16::MAX)
            .to_le_bytes(),
    );
    out.extend_from_slice(&compressed);

    let mut hasher = Hasher::new();
    hasher.update(input);
    out.extend_from_slice(&hasher.finalize().to_le_bytes());
    out.extend_from_slice(&(input.len() as u32).to_le_bytes());
    Ok(out)
}

#[cfg(test)]
mod tests {
    use std::io::{Read, Write};

    use super::*;

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
}
