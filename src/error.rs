use std::fmt;

pub type Result<T> = std::result::Result<T, FastqError>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqPosition {
    pub byte_offset: u64,
    pub record_index: u64,
    pub line_index: u8,
}

impl FastqPosition {
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

#[derive(Debug)]
#[non_exhaustive]
pub enum FastqError {
    Io(std::io::Error),
    Format(String),
    FormatAt {
        message: String,
        position: FastqPosition,
    },
    Bgzf(String),
    RecordTooLarge {
        slab_size: usize,
    },
}

impl fmt::Display for FastqError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(err) => write!(f, "I/O error: {err}"),
            Self::Format(msg) => write!(f, "FASTQ parse error: {msg}"),
            Self::FormatAt { message, position } => {
                write!(f, "FASTQ parse error at {position}: {message}")
            }
            Self::Bgzf(msg) => write!(f, "BGZF error: {msg}"),
            Self::RecordTooLarge { slab_size } => {
                write!(f, "FASTQ record exceeds slab size ({slab_size} bytes)")
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
