use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{Receiver, SyncSender, TrySendError, sync_channel};
use std::thread::{self, JoinHandle};
use std::time::Duration;

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

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BgzfVirtualOffset(u64);

impl BgzfVirtualOffset {
    const MAX_COMPRESSED_OFFSET: u64 = (1_u64 << 48) - 1;

    pub fn from_parts(compressed_offset: u64, in_block_offset: u16) -> Result<Self> {
        if compressed_offset > Self::MAX_COMPRESSED_OFFSET {
            return Err(FastqError::Bgzf(
                "BGZF compressed offset exceeds virtual-offset range".into(),
            ));
        }
        Ok(Self((compressed_offset << 16) | u64::from(in_block_offset)))
    }

    pub fn from_raw(raw: u64) -> Self {
        Self(raw)
    }

    pub fn raw(self) -> u64 {
        self.0
    }

    pub fn compressed_offset(self) -> u64 {
        self.0 >> 16
    }

    pub fn in_block_offset(self) -> u16 {
        (self.0 & 0xffff) as u16
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BgzfIndexEntry {
    pub compressed_offset: u64,
    pub uncompressed_offset: u64,
    pub compressed_size: u32,
    pub uncompressed_size: u32,
}

impl BgzfIndexEntry {
    pub fn block_virtual_offset(&self) -> Result<BgzfVirtualOffset> {
        BgzfVirtualOffset::from_parts(self.compressed_offset, 0)
    }

    pub fn contains_uncompressed_offset(&self, offset: u64) -> bool {
        offset >= self.uncompressed_offset
            && offset < self.uncompressed_offset + u64::from(self.uncompressed_size)
    }

    pub fn virtual_offset_for(&self, offset: u64) -> Result<Option<BgzfVirtualOffset>> {
        if !self.contains_uncompressed_offset(offset) {
            return Ok(None);
        }
        let in_block = offset - self.uncompressed_offset;
        let in_block = u16::try_from(in_block)
            .map_err(|_| FastqError::Bgzf("BGZF in-block offset exceeds u16 range".into()))?;
        BgzfVirtualOffset::from_parts(self.compressed_offset, in_block).map(Some)
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BgzfIndex {
    entries: Vec<BgzfIndexEntry>,
    uncompressed_len: u64,
    compressed_len: u64,
}

impl BgzfIndex {
    pub fn entries(&self) -> &[BgzfIndexEntry] {
        &self.entries
    }

    pub fn uncompressed_len(&self) -> u64 {
        self.uncompressed_len
    }

    pub fn compressed_len(&self) -> u64 {
        self.compressed_len
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn entry_for_uncompressed_offset(&self, offset: u64) -> Option<&BgzfIndexEntry> {
        if offset >= self.uncompressed_len {
            return None;
        }
        let idx = self
            .entries
            .partition_point(|entry| entry.uncompressed_offset <= offset);
        idx.checked_sub(1)
            .and_then(|entry_idx| self.entries.get(entry_idx))
            .filter(|entry| entry.contains_uncompressed_offset(offset))
    }

    pub fn virtual_offset_for_uncompressed_offset(
        &self,
        offset: u64,
    ) -> Result<Option<BgzfVirtualOffset>> {
        let Some(entry) = self.entry_for_uncompressed_offset(offset) else {
            return Ok(None);
        };
        entry.virtual_offset_for(offset)
    }
}

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

pub struct BgzfParallelReader {
    result_rx: Receiver<ParallelMsg>,
    current: Vec<u8>,
    pos: usize,
    pending: BTreeMap<usize, Vec<u8>>,
    next_index: usize,
    total_blocks: Option<usize>,
    finished: bool,
    cancel: Arc<AtomicBool>,
    handles: Vec<JoinHandle<()>>,
}

enum Job {
    Block(usize, CompressedBlock),
    End,
}

enum ParallelMsg {
    Data(usize, Result<Vec<u8>>),
    End(usize),
    Fatal(String),
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

impl BgzfParallelReader {
    pub fn new<R>(inner: R, workers: usize) -> Result<Self>
    where
        R: Read + Send + 'static,
    {
        let worker_count = workers.max(1);
        let cancel = Arc::new(AtomicBool::new(false));
        let (result_tx, result_rx) = sync_channel(worker_count.saturating_mul(2).max(2));

        let mut job_txs = Vec::with_capacity(worker_count);
        let mut job_rxs = Vec::with_capacity(worker_count);
        for _ in 0..worker_count {
            let (tx, rx) = sync_channel(2);
            job_txs.push(tx);
            job_rxs.push(rx);
        }

        let mut handles = Vec::with_capacity(worker_count + 1);
        for (worker_id, rx) in job_rxs.into_iter().enumerate() {
            let tx = result_tx.clone();
            let cancel_worker = Arc::clone(&cancel);
            let handle = thread::Builder::new()
                .name(format!("microraptor-bgzf-decode-{worker_id}"))
                .spawn(move || bgzf_worker_loop(rx, tx, cancel_worker))?;
            handles.push(handle);
        }

        let cancel_reader = Arc::clone(&cancel);
        let handle = thread::Builder::new()
            .name("microraptor-bgzf-read".into())
            .spawn(move || bgzf_reader_loop(inner, job_txs, result_tx, cancel_reader))?;
        handles.push(handle);

        Ok(Self {
            result_rx,
            current: Vec::new(),
            pos: 0,
            pending: BTreeMap::new(),
            next_index: 0,
            total_blocks: None,
            finished: false,
            cancel,
            handles,
        })
    }

    fn refill_ordered(&mut self) -> std::io::Result<bool> {
        loop {
            if let Some(buf) = self.pending.remove(&self.next_index) {
                self.current = buf;
                self.pos = 0;
                self.next_index += 1;
                if !self.current.is_empty() {
                    return Ok(true);
                }
                continue;
            }

            if self.total_blocks == Some(self.next_index) {
                self.finished = true;
                return Ok(false);
            }

            match self.result_rx.recv() {
                Ok(ParallelMsg::Data(index, Ok(buf))) => {
                    if index == self.next_index {
                        self.current = buf;
                        self.pos = 0;
                        self.next_index += 1;
                        if !self.current.is_empty() {
                            return Ok(true);
                        }
                    } else {
                        self.pending.insert(index, buf);
                    }
                }
                Ok(ParallelMsg::Data(_, Err(err))) => {
                    self.finished = true;
                    return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, err));
                }
                Ok(ParallelMsg::End(total)) => {
                    self.total_blocks = Some(total);
                }
                Ok(ParallelMsg::Fatal(msg)) => {
                    self.finished = true;
                    return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, msg));
                }
                Err(_) => {
                    self.finished = true;
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::UnexpectedEof,
                        "BGZF parallel pipeline ended before EOF",
                    ));
                }
            }
        }
    }
}

impl Read for BgzfParallelReader {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        if self.pos >= self.current.len() {
            self.current.clear();
            self.pos = 0;
            if self.finished || !self.refill_ordered()? {
                return Ok(0);
            }
        }
        let n = out.len().min(self.current.len() - self.pos);
        out[..n].copy_from_slice(&self.current[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

impl Drop for BgzfParallelReader {
    fn drop(&mut self) {
        self.cancel.store(true, Ordering::Release);
        for handle in self.handles.drain(..) {
            let _ = handle.join();
        }
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

pub fn build_bgzf_index<R: Read>(mut reader: R) -> Result<BgzfIndex> {
    let mut entries = Vec::new();
    let mut compressed_offset = 0_u64;
    let mut uncompressed_offset = 0_u64;

    while let Some(block) = read_block(&mut reader)? {
        let compressed_size = u32::try_from(block.bytes.len())
            .map_err(|_| FastqError::Bgzf("BGZF block size exceeds u32 range".into()))?;
        if block.is_eof() {
            compressed_offset += u64::from(compressed_size);
            break;
        }

        let decoded = decode_block(&block)?;
        let uncompressed_size = u32::try_from(decoded.len()).map_err(|_| {
            FastqError::Bgzf("BGZF uncompressed block size exceeds u32 range".into())
        })?;
        if uncompressed_size > 65_536 {
            return Err(FastqError::Bgzf(
                "BGZF uncompressed block exceeds 64 KiB".into(),
            ));
        }

        entries.push(BgzfIndexEntry {
            compressed_offset,
            uncompressed_offset,
            compressed_size,
            uncompressed_size,
        });
        compressed_offset += u64::from(compressed_size);
        uncompressed_offset += u64::from(uncompressed_size);
    }

    Ok(BgzfIndex {
        entries,
        uncompressed_len: uncompressed_offset,
        compressed_len: compressed_offset,
    })
}

fn bgzf_reader_loop<R>(
    mut inner: R,
    job_txs: Vec<SyncSender<Job>>,
    result_tx: SyncSender<ParallelMsg>,
    cancel: Arc<AtomicBool>,
) where
    R: Read,
{
    let mut index = 0;
    let mut fatal = None;
    while !cancel.load(Ordering::Acquire) {
        match read_block(&mut inner) {
            Ok(Some(block)) if block.is_eof() => break,
            Ok(Some(block)) => {
                let worker = index % job_txs.len();
                if !send_job(&job_txs[worker], Job::Block(index, block), &cancel) {
                    return;
                }
                index += 1;
            }
            Ok(None) => break,
            Err(err) => {
                fatal = Some(err.to_string());
                break;
            }
        }
    }

    for tx in &job_txs {
        let _ = send_job(tx, Job::End, &cancel);
    }

    if let Some(msg) = fatal {
        let _ = send_result(&result_tx, ParallelMsg::Fatal(msg), &cancel);
    } else {
        let _ = send_result(&result_tx, ParallelMsg::End(index), &cancel);
    }
}

fn bgzf_worker_loop(
    rx: Receiver<Job>,
    result_tx: SyncSender<ParallelMsg>,
    cancel: Arc<AtomicBool>,
) {
    while !cancel.load(Ordering::Acquire) {
        match rx.recv_timeout(Duration::from_millis(10)) {
            Ok(Job::Block(index, block)) => {
                let decoded = decode_block(&block);
                if !send_result(&result_tx, ParallelMsg::Data(index, decoded), &cancel) {
                    return;
                }
            }
            Ok(Job::End) => return,
            Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {}
            Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => return,
        }
    }
}

fn send_job(tx: &SyncSender<Job>, msg: Job, cancel: &AtomicBool) -> bool {
    send_bounded(tx, msg, cancel)
}

fn send_result(tx: &SyncSender<ParallelMsg>, msg: ParallelMsg, cancel: &AtomicBool) -> bool {
    send_bounded(tx, msg, cancel)
}

fn send_bounded<T>(tx: &SyncSender<T>, mut msg: T, cancel: &AtomicBool) -> bool {
    while !cancel.load(Ordering::Acquire) {
        match tx.try_send(msg) {
            Ok(()) => return true,
            Err(TrySendError::Full(returned)) => {
                msg = returned;
                thread::sleep(Duration::from_millis(1));
            }
            Err(TrySendError::Disconnected(_)) => return false,
        }
    }
    false
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
    fn virtual_offset_checks_compressed_offset_range() {
        let err = BgzfVirtualOffset::from_parts(1_u64 << 48, 0).unwrap_err();
        assert!(err.to_string().contains("virtual-offset range"));

        let vo = BgzfVirtualOffset::from_parts(123, 45).unwrap();
        assert_eq!(vo.raw(), (123 << 16) | 45);
        assert_eq!(vo.compressed_offset(), 123);
        assert_eq!(vo.in_block_offset(), 45);
    }
}
