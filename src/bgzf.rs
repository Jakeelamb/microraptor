//! BGZF readers, writers, indexing, and parallel block helpers.
//!
//! BGZF is a blocked gzip variant used by htslib-compatible bioinformatics
//! tools. Microraptor exposes these types so FASTQ callers can use the same
//! parser over raw, ordinary gzip, and BGZF transports.

use std::collections::BTreeMap;
use std::io::{Read, Seek, SeekFrom, Write};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
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
const DEFAULT_PARALLEL_MIN_COMPRESSED_BYTES: u64 = 512 * 1024;

/// Canonical empty BGZF EOF marker block.
pub const BGZF_EOF_BLOCK: &[u8] = &[
    31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 27, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

/// Inflate backend used for BGZF block decompression.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BgzfInflateBackend {
    /// Use the `flate2` backend.
    Flate2,
    /// Use libdeflate when the `libdeflate` feature is enabled.
    #[cfg(feature = "libdeflate")]
    Libdeflate,
}

impl Default for BgzfInflateBackend {
    fn default() -> Self {
        Self::fastest_available()
    }
}

/// Deflate backend used for BGZF block compression.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BgzfDeflateBackend {
    /// Use the `flate2` backend.
    Flate2,
    /// Use libdeflate when the `libdeflate` feature is enabled.
    #[cfg(feature = "libdeflate")]
    Libdeflate,
}

impl Default for BgzfDeflateBackend {
    fn default() -> Self {
        Self::fastest_available()
    }
}

impl BgzfDeflateBackend {
    /// Return the fastest compiled-in deflate backend.
    pub const fn fastest_available() -> Self {
        #[cfg(feature = "libdeflate")]
        {
            Self::Libdeflate
        }
        #[cfg(not(feature = "libdeflate"))]
        {
            Self::Flate2
        }
    }
}

/// Configuration for bounded parallel BGZF decoding.
#[derive(Debug, Clone)]
pub struct BgzfParallelConfig {
    /// Number of decode workers.
    pub workers: usize,
    /// Bounded queue depth from reader to worker threads.
    pub job_queue_depth: usize,
    /// Bounded queue depth from worker threads to the ordered output reader.
    pub result_queue_depth: usize,
    /// Inflate backend used by workers.
    pub backend: BgzfInflateBackend,
    /// Minimum compressed input size before adaptive readers choose parallel mode.
    pub parallel_min_compressed_bytes: u64,
    /// Optional backpressure metrics collector.
    pub metrics: Option<Arc<BgzfPipelineMetrics>>,
}

impl Default for BgzfParallelConfig {
    fn default() -> Self {
        Self {
            workers: 1,
            job_queue_depth: 2,
            result_queue_depth: 2,
            backend: BgzfInflateBackend::default(),
            parallel_min_compressed_bytes: DEFAULT_PARALLEL_MIN_COMPRESSED_BYTES,
            metrics: None,
        }
    }
}

/// Shared counters for BGZF pipeline backpressure observations.
#[derive(Debug, Default)]
pub struct BgzfPipelineMetrics {
    job_queue_full: AtomicU64,
    result_queue_full: AtomicU64,
}

/// Snapshot of BGZF pipeline backpressure counters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BgzfPipelineMetricsSnapshot {
    /// Number of times the reader observed a full job queue.
    pub job_queue_full: u64,
    /// Number of times a worker observed a full result queue.
    pub result_queue_full: u64,
}

impl BgzfPipelineMetrics {
    /// Read the current metric counters.
    pub fn snapshot(&self) -> BgzfPipelineMetricsSnapshot {
        BgzfPipelineMetricsSnapshot {
            job_queue_full: self.job_queue_full.load(Ordering::Relaxed),
            result_queue_full: self.result_queue_full.load(Ordering::Relaxed),
        }
    }

    fn record(&self, channel: BgzfBackpressureChannel) {
        match channel {
            BgzfBackpressureChannel::Job => {
                self.job_queue_full.fetch_add(1, Ordering::Relaxed);
            }
            BgzfBackpressureChannel::Result => {
                self.result_queue_full.fetch_add(1, Ordering::Relaxed);
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum BgzfBackpressureChannel {
    Job,
    Result,
}

impl BgzfParallelConfig {
    /// Construct a parallel config for `workers` decode workers.
    pub fn new(workers: usize) -> Self {
        Self {
            workers,
            result_queue_depth: workers.saturating_mul(4).max(2),
            ..Self::default()
        }
    }

    /// Set the inflate backend.
    pub fn with_inflate_backend(mut self, backend: BgzfInflateBackend) -> Self {
        self.backend = backend;
        self
    }

    /// Set bounded queue depths.
    pub fn with_queue_depths(mut self, job_queue_depth: usize, result_queue_depth: usize) -> Self {
        self.job_queue_depth = job_queue_depth.max(1);
        self.result_queue_depth = result_queue_depth.max(1);
        self
    }

    /// Set the adaptive serial/parallel compressed-size threshold.
    pub fn with_parallel_min_compressed_bytes(mut self, bytes: u64) -> Self {
        self.parallel_min_compressed_bytes = bytes;
        self
    }

    /// Attach a shared metrics collector.
    pub fn with_metrics(mut self, metrics: Arc<BgzfPipelineMetrics>) -> Self {
        self.metrics = Some(metrics);
        self
    }

    /// Return true when this config should parallelize an input of `compressed_len`.
    pub fn should_parallelize(&self, compressed_len: u64) -> bool {
        self.workers > 1 && compressed_len >= self.parallel_min_compressed_bytes
    }
}

impl BgzfInflateBackend {
    /// Return the fastest compiled-in inflate backend.
    pub const fn fastest_available() -> Self {
        #[cfg(feature = "libdeflate")]
        {
            Self::Libdeflate
        }
        #[cfg(not(feature = "libdeflate"))]
        {
            Self::Flate2
        }
    }
}

/// BGZF virtual offset.
///
/// The upper 48 bits are the compressed block offset, and the lower 16 bits are
/// the uncompressed offset inside that block.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BgzfVirtualOffset(u64);

impl BgzfVirtualOffset {
    const MAX_COMPRESSED_OFFSET: u64 = (1_u64 << 48) - 1;

    /// Build a virtual offset from compressed and in-block offsets.
    pub fn from_parts(compressed_offset: u64, in_block_offset: u16) -> Result<Self> {
        if compressed_offset > Self::MAX_COMPRESSED_OFFSET {
            return Err(FastqError::Bgzf(
                "BGZF compressed offset exceeds virtual-offset range".into(),
            ));
        }
        Ok(Self((compressed_offset << 16) | u64::from(in_block_offset)))
    }

    /// Wrap a raw virtual-offset integer.
    pub fn from_raw(raw: u64) -> Self {
        Self(raw)
    }

    /// Return the raw virtual-offset integer.
    pub fn raw(self) -> u64 {
        self.0
    }

    /// Return the compressed block byte offset.
    pub fn compressed_offset(self) -> u64 {
        self.0 >> 16
    }

    /// Return the uncompressed offset inside the BGZF block.
    pub fn in_block_offset(self) -> u16 {
        (self.0 & 0xffff) as u16
    }
}

/// One BGZF index entry for a compressed block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BgzfIndexEntry {
    /// Compressed stream byte offset of the block.
    pub compressed_offset: u64,
    /// Uncompressed stream byte offset of the block.
    pub uncompressed_offset: u64,
    /// Compressed block size in bytes.
    pub compressed_size: u32,
    /// Uncompressed block size in bytes.
    pub uncompressed_size: u32,
}

impl BgzfIndexEntry {
    /// Return the virtual offset at the start of this block.
    pub fn block_virtual_offset(&self) -> Result<BgzfVirtualOffset> {
        BgzfVirtualOffset::from_parts(self.compressed_offset, 0)
    }

    /// Return true if an uncompressed stream offset falls inside this block.
    pub fn contains_uncompressed_offset(&self, offset: u64) -> bool {
        offset >= self.uncompressed_offset
            && offset < self.uncompressed_offset + u64::from(self.uncompressed_size)
    }

    /// Convert an uncompressed stream offset to a virtual offset if it is inside this block.
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

/// In-memory index of BGZF block offsets.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BgzfIndex {
    entries: Vec<BgzfIndexEntry>,
    uncompressed_len: u64,
    compressed_len: u64,
}

impl BgzfIndex {
    /// Return all block index entries.
    pub fn entries(&self) -> &[BgzfIndexEntry] {
        &self.entries
    }

    /// Return total uncompressed stream length represented by the index.
    pub fn uncompressed_len(&self) -> u64 {
        self.uncompressed_len
    }

    /// Return total compressed stream length consumed while building the index.
    pub fn compressed_len(&self) -> u64 {
        self.compressed_len
    }

    /// Return true when the index has no data blocks.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Return the number of indexed data blocks.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Find the block containing an uncompressed stream offset.
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

    /// Convert an uncompressed stream offset into a BGZF virtual offset.
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

/// Serial BGZF reader implementing [`std::io::Read`].
pub struct BgzfReader<R> {
    inner: R,
    backend: BgzfInflateBackend,
    decoded: Vec<u8>,
    pos: usize,
    eof: bool,
}

/// Seekable BGZF reader using virtual offsets.
pub struct BgzfSeekReader<R> {
    inner: R,
    backend: BgzfInflateBackend,
    decoded: Vec<u8>,
    pos: usize,
    eof: bool,
}

/// Bounded parallel BGZF reader preserving block order.
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

/// Adaptive BGZF reader that selects serial or parallel decoding.
pub enum BgzfAutoReader<R> {
    /// Serial reader variant.
    Serial(BgzfReader<R>),
    /// Parallel reader variant.
    Parallel(BgzfParallelReader),
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
    /// Construct a serial BGZF reader using the default inflate backend.
    pub fn new(inner: R) -> Self {
        Self::with_inflate_backend(inner, BgzfInflateBackend::default())
    }

    /// Construct a serial BGZF reader using an explicit inflate backend.
    pub fn with_inflate_backend(inner: R, backend: BgzfInflateBackend) -> Self {
        Self {
            inner,
            backend,
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
            decode_block_into_with_backend(&block, self.backend, &mut self.decoded)
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

impl<R: Read + Seek> BgzfSeekReader<R> {
    /// Construct a seekable BGZF reader using the default inflate backend.
    pub fn new(inner: R) -> Self {
        Self::with_inflate_backend(inner, BgzfInflateBackend::default())
    }

    /// Construct a seekable BGZF reader using an explicit inflate backend.
    pub fn with_inflate_backend(inner: R, backend: BgzfInflateBackend) -> Self {
        Self {
            inner,
            backend,
            decoded: Vec::new(),
            pos: 0,
            eof: false,
        }
    }

    /// Seek to a BGZF virtual offset.
    pub fn seek_virtual_offset(&mut self, offset: BgzfVirtualOffset) -> std::io::Result<()> {
        self.inner
            .seek(SeekFrom::Start(offset.compressed_offset()))?;
        self.decoded.clear();
        self.pos = 0;
        self.eof = false;

        let Some(block) = read_block(&mut self.inner)? else {
            return Err(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "BGZF virtual offset points past end of stream",
            ));
        };

        if block.is_eof() {
            if offset.in_block_offset() == 0 {
                self.eof = true;
                return Ok(());
            }
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "BGZF virtual offset points inside EOF block",
            ));
        }

        decode_block_into_with_backend(&block, self.backend, &mut self.decoded)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        let pos = usize::from(offset.in_block_offset());
        if pos > self.decoded.len() {
            self.decoded.clear();
            self.pos = 0;
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "BGZF in-block virtual offset exceeds decoded block length",
            ));
        }
        self.pos = pos;
        Ok(())
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
            decode_block_into_with_backend(&block, self.backend, &mut self.decoded)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            if !self.decoded.is_empty() {
                return Ok(());
            }
        }
    }
}

impl<R: Read + Seek> Read for BgzfSeekReader<R> {
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
    /// Construct a parallel BGZF reader with `workers` decode workers.
    pub fn new<R>(inner: R, workers: usize) -> Result<Self>
    where
        R: Read + Send + 'static,
    {
        Self::with_inflate_backend(inner, workers, BgzfInflateBackend::default())
    }

    /// Construct a parallel BGZF reader with an explicit inflate backend.
    pub fn with_inflate_backend<R>(
        inner: R,
        workers: usize,
        backend: BgzfInflateBackend,
    ) -> Result<Self>
    where
        R: Read + Send + 'static,
    {
        Self::with_config(
            inner,
            BgzfParallelConfig::new(workers).with_inflate_backend(backend),
        )
    }

    /// Construct a parallel BGZF reader with full pipeline configuration.
    pub fn with_config<R>(inner: R, config: BgzfParallelConfig) -> Result<Self>
    where
        R: Read + Send + 'static,
    {
        let worker_count = config.workers.max(1);
        let job_queue_depth = config.job_queue_depth.max(1);
        let result_queue_depth = config.result_queue_depth.max(1);
        let backend = config.backend;
        let metrics = config.metrics.clone();
        let cancel = Arc::new(AtomicBool::new(false));
        let (result_tx, result_rx) = sync_channel(result_queue_depth);

        let mut job_txs = Vec::with_capacity(worker_count);
        let mut job_rxs = Vec::with_capacity(worker_count);
        for _ in 0..worker_count {
            let (tx, rx) = sync_channel(job_queue_depth);
            job_txs.push(tx);
            job_rxs.push(rx);
        }

        let mut handles = Vec::with_capacity(worker_count + 1);
        for (worker_id, rx) in job_rxs.into_iter().enumerate() {
            let tx = result_tx.clone();
            let cancel_worker = Arc::clone(&cancel);
            let metrics_worker = metrics.clone();
            let handle = thread::Builder::new()
                .name(format!("microraptor-bgzf-decode-{worker_id}"))
                .spawn(move || bgzf_worker_loop(rx, tx, cancel_worker, backend, metrics_worker))?;
            handles.push(handle);
        }

        let cancel_reader = Arc::clone(&cancel);
        let metrics_reader = metrics.clone();
        let handle = thread::Builder::new()
            .name("microraptor-bgzf-read".into())
            .spawn(move || {
                bgzf_reader_loop(inner, job_txs, result_tx, cancel_reader, metrics_reader)
            })?;
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

impl<R> BgzfAutoReader<R>
where
    R: Read + Send + 'static,
{
    /// Construct an adaptive BGZF reader from compressed input length and config.
    pub fn with_config(inner: R, compressed_len: u64, config: BgzfParallelConfig) -> Result<Self> {
        if config.should_parallelize(compressed_len) {
            Ok(Self::Parallel(BgzfParallelReader::with_config(
                inner, config,
            )?))
        } else {
            Ok(Self::Serial(BgzfReader::with_inflate_backend(
                inner,
                config.backend,
            )))
        }
    }
}

impl<R: Read> Read for BgzfAutoReader<R> {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        match self {
            Self::Serial(reader) => reader.read(out),
            Self::Parallel(reader) => reader.read(out),
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

/// BGZF writer implementing [`std::io::Write`].
pub struct BgzfWriter<W> {
    inner: Option<W>,
    pending: Vec<u8>,
    level: Compression,
    backend: BgzfDeflateBackend,
}

impl<W: Write> BgzfWriter<W> {
    /// Construct a BGZF writer using fast compression.
    pub fn new(inner: W) -> Self {
        Self::with_compression(inner, Compression::fast())
    }

    /// Construct a BGZF writer using an explicit flate2 compression level.
    pub fn with_compression(inner: W, level: Compression) -> Self {
        Self {
            inner: Some(inner),
            pending: Vec::with_capacity(BGZF_MAX_PAYLOAD),
            level,
            backend: BgzfDeflateBackend::Flate2,
        }
    }

    /// Construct a BGZF writer using an explicit deflate backend.
    pub fn with_deflate_backend(inner: W, backend: BgzfDeflateBackend) -> Self {
        Self {
            inner: Some(inner),
            pending: Vec::with_capacity(BGZF_MAX_PAYLOAD),
            level: Compression::fast(),
            backend,
        }
    }

    /// Finish the stream, write the EOF block, flush, and return the inner writer.
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
        let block = encode_block_with_backend(&self.pending, self.level, self.backend)?;
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

/// Decompress a complete BGZF stream using parallel block decoding.
pub fn decompress_bgzf_parallel<R: Read>(reader: R, workers: usize) -> Result<Vec<u8>> {
    decompress_bgzf_parallel_with_inflate_backend(reader, workers, BgzfInflateBackend::default())
}

/// Decompress a complete BGZF stream using an explicit inflate backend.
pub fn decompress_bgzf_parallel_with_inflate_backend<R: Read>(
    mut reader: R,
    workers: usize,
    backend: BgzfInflateBackend,
) -> Result<Vec<u8>> {
    let mut blocks = Vec::new();
    while let Some(block) = read_block(&mut reader)? {
        if block.is_eof() {
            break;
        }
        blocks.push(block);
    }
    let decoded = parallel_map(blocks, workers, |block| {
        decode_block_with_backend(&block, backend)
    })?;
    let total = decoded.iter().map(Vec::len).sum();
    let mut out = Vec::with_capacity(total);
    for chunk in decoded {
        out.extend_from_slice(&chunk);
    }
    Ok(out)
}

/// Compress a complete buffer as BGZF using parallel block compression.
pub fn compress_bgzf_parallel(input: &[u8], workers: usize) -> Result<Vec<u8>> {
    compress_bgzf_parallel_with_deflate_backend(input, workers, BgzfDeflateBackend::Flate2)
}

/// Compress a complete buffer as BGZF using an explicit deflate backend.
pub fn compress_bgzf_parallel_with_deflate_backend(
    input: &[u8],
    workers: usize,
    backend: BgzfDeflateBackend,
) -> Result<Vec<u8>> {
    let chunks: Vec<&[u8]> = input.chunks(BGZF_MAX_PAYLOAD).collect();
    let encoded = parallel_map(chunks, workers, |chunk| {
        encode_block_with_backend(chunk, Compression::fast(), backend)
    })?;
    let total = encoded.iter().map(Vec::len).sum::<usize>() + BGZF_EOF_BLOCK.len();
    let mut out = Vec::with_capacity(total);
    for block in encoded {
        out.extend_from_slice(&block);
    }
    out.extend_from_slice(BGZF_EOF_BLOCK);
    Ok(out)
}

/// Build a BGZF block index from a complete BGZF stream.
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
    metrics: Option<Arc<BgzfPipelineMetrics>>,
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
                if !send_job(
                    &job_txs[worker],
                    Job::Block(index, block),
                    &cancel,
                    metrics.as_deref(),
                ) {
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
        let _ = send_job(tx, Job::End, &cancel, metrics.as_deref());
    }

    if let Some(msg) = fatal {
        let _ = send_result(
            &result_tx,
            ParallelMsg::Fatal(msg),
            &cancel,
            metrics.as_deref(),
        );
    } else {
        let _ = send_result(
            &result_tx,
            ParallelMsg::End(index),
            &cancel,
            metrics.as_deref(),
        );
    }
}

fn bgzf_worker_loop(
    rx: Receiver<Job>,
    result_tx: SyncSender<ParallelMsg>,
    cancel: Arc<AtomicBool>,
    backend: BgzfInflateBackend,
    metrics: Option<Arc<BgzfPipelineMetrics>>,
) {
    while !cancel.load(Ordering::Acquire) {
        match rx.recv_timeout(Duration::from_millis(10)) {
            Ok(Job::Block(index, block)) => {
                let decoded = decode_block_with_backend(&block, backend);
                if !send_result(
                    &result_tx,
                    ParallelMsg::Data(index, decoded),
                    &cancel,
                    metrics.as_deref(),
                ) {
                    return;
                }
            }
            Ok(Job::End) => return,
            Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {}
            Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => return,
        }
    }
}

fn send_job(
    tx: &SyncSender<Job>,
    msg: Job,
    cancel: &AtomicBool,
    metrics: Option<&BgzfPipelineMetrics>,
) -> bool {
    send_bounded(tx, msg, cancel, metrics, BgzfBackpressureChannel::Job)
}

fn send_result(
    tx: &SyncSender<ParallelMsg>,
    msg: ParallelMsg,
    cancel: &AtomicBool,
    metrics: Option<&BgzfPipelineMetrics>,
) -> bool {
    send_bounded(tx, msg, cancel, metrics, BgzfBackpressureChannel::Result)
}

fn send_bounded<T>(
    tx: &SyncSender<T>,
    mut msg: T,
    cancel: &AtomicBool,
    metrics: Option<&BgzfPipelineMetrics>,
    channel: BgzfBackpressureChannel,
) -> bool {
    while !cancel.load(Ordering::Acquire) {
        match tx.try_send(msg) {
            Ok(()) => return true,
            Err(TrySendError::Full(returned)) => {
                if let Some(metrics) = metrics {
                    metrics.record(channel);
                }
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
    match reader.read(&mut header[..1]) {
        Ok(0) => return Ok(None),
        Ok(_) => {}
        Err(e) => return Err(e),
    }
    reader.read_exact(&mut header[1..])?;
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
    decode_block_with_backend(block, BgzfInflateBackend::default())
}

fn decode_block_with_backend(
    block: &CompressedBlock,
    backend: BgzfInflateBackend,
) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    decode_block_into_with_backend(block, backend, &mut out)?;
    Ok(out)
}

fn decode_block_into_with_backend(
    block: &CompressedBlock,
    backend: BgzfInflateBackend,
    out: &mut Vec<u8>,
) -> Result<()> {
    let bytes = &block.bytes;
    if bytes.len() < BGZF_HEADER_LEN + GZIP_TRAILER_LEN {
        return Err(FastqError::Bgzf("short block".into()));
    }
    let compressed_end = bytes.len() - GZIP_TRAILER_LEN;

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

    let deflate = &bytes[BGZF_HEADER_LEN..compressed_end];
    match backend {
        BgzfInflateBackend::Flate2 => inflate_block_flate2_into(deflate, expected_len, out)?,
        #[cfg(feature = "libdeflate")]
        BgzfInflateBackend::Libdeflate => {
            inflate_block_libdeflate_into(deflate, expected_len, out)?
        }
    }
    if out.len() != expected_len {
        return Err(FastqError::Bgzf("uncompressed size mismatch".into()));
    }
    let mut hasher = Hasher::new();
    hasher.update(out);
    if hasher.finalize() != expected_crc {
        return Err(FastqError::Bgzf("CRC32 mismatch".into()));
    }
    Ok(())
}

fn inflate_block_flate2_into(deflate: &[u8], expected_len: usize, out: &mut Vec<u8>) -> Result<()> {
    out.clear();
    if expected_len > 0 {
        out.reserve(expected_len.saturating_sub(out.capacity()));
    }
    let mut decoder = DeflateDecoder::new(deflate);
    decoder.read_to_end(out)?;
    Ok(())
}

#[cfg(feature = "libdeflate")]
fn inflate_block_libdeflate_into(
    deflate: &[u8],
    expected_len: usize,
    out: &mut Vec<u8>,
) -> Result<()> {
    out.resize(expected_len, 0);
    let mut decompressor = libdeflater::Decompressor::new();
    let actual = decompressor
        .deflate_decompress(deflate, out)
        .map_err(|e| FastqError::Bgzf(format!("libdeflate inflate failed: {e}")))?;
    out.truncate(actual);
    Ok(())
}

fn encode_block_with_backend(
    input: &[u8],
    level: Compression,
    backend: BgzfDeflateBackend,
) -> Result<Vec<u8>> {
    if input.len() > BGZF_MAX_PAYLOAD {
        return Err(FastqError::Bgzf("payload exceeds BGZF block policy".into()));
    }

    let compressed = match backend {
        BgzfDeflateBackend::Flate2 => deflate_block_flate2(input, level)?,
        #[cfg(feature = "libdeflate")]
        BgzfDeflateBackend::Libdeflate => deflate_block_libdeflate(input)?,
    };
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

fn deflate_block_flate2(input: &[u8], level: Compression) -> Result<Vec<u8>> {
    let mut encoder = DeflateEncoder::new(Vec::new(), level);
    encoder.write_all(input)?;
    Ok(encoder.finish()?)
}

#[cfg(feature = "libdeflate")]
fn deflate_block_libdeflate(input: &[u8]) -> Result<Vec<u8>> {
    let mut compressor = libdeflater::Compressor::new(libdeflater::CompressionLvl::fastest());
    let mut out = vec![0_u8; compressor.deflate_compress_bound(input.len())];
    let n = compressor
        .deflate_compress(input, &mut out)
        .map_err(|e| FastqError::Bgzf(format!("libdeflate deflate failed: {e}")))?;
    out.truncate(n);
    Ok(out)
}

#[cfg(test)]
mod tests;
