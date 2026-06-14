use std::fmt;

/// Crate-local result type.
pub type Result<T> = std::result::Result<T, FastqError>;

/// Location of a FASTQ or FASTA parsing error.
///
/// `byte_offset` is absolute within the original byte stream. `record_index`
/// is zero-based. `line_index` is `0` for header, `1` for sequence, `2` for
/// plus, and `3` for quality.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqPosition {
    /// Absolute byte offset in the original input stream.
    pub byte_offset: u64,
    /// Zero-based FASTQ record index.
    pub record_index: u64,
    /// Zero-based line within the four-line FASTQ record.
    pub line_index: u8,
}

impl FastqPosition {
    /// Create a new position.
    pub const fn new(byte_offset: u64, record_index: u64, line_index: u8) -> Self {
        Self {
            byte_offset,
            record_index,
            line_index,
        }
    }
}

impl fmt::Display for FastqPosition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "byte {}, record {}, line {}",
            self.byte_offset, self.record_index, self.line_index
        )
    }
}

/// Error type for FASTQ, FASTA, gzip, and BGZF operations.
#[derive(Debug)]
#[non_exhaustive]
pub enum FastqError {
    /// I/O error from the underlying reader or file.
    Io(std::io::Error),
    /// Format error without a precise byte position.
    Format(String),
    /// Format error with byte, record, and line position.
    FormatAt {
        /// Human-readable format error.
        message: String,
        /// Position of the error.
        position: FastqPosition,
    },
    /// BGZF-specific error.
    Bgzf(String),
    /// A single record did not fit in the configured slab.
    RecordTooLarge {
        /// Configured slab size in bytes.
        slab_size: usize,
    },
}

impl fmt::Display for FastqError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(err) => write!(f, "I/O error: {err}"),
            Self::Format(msg) => write!(f, "parse error: {msg}"),
            Self::FormatAt { message, position } => {
                write!(f, "parse error at {position}: {message}")
            }
            Self::Bgzf(msg) => write!(f, "BGZF error: {msg}"),
            Self::RecordTooLarge { slab_size } => {
                write!(f, "record exceeds slab size ({slab_size} bytes)")
            }
        }
    }
}

impl std::error::Error for FastqError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(err) => Some(err),
            Self::Format(_)
            | Self::FormatAt { .. }
            | Self::Bgzf(_)
            | Self::RecordTooLarge { .. } => None,
        }
    }
}

impl From<std::io::Error> for FastqError {
    fn from(value: std::io::Error) -> Self {
        Self::Io(value)
    }
}
