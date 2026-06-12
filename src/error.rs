use std::fmt;

pub type Result<T> = std::result::Result<T, FastqError>;

#[derive(Debug)]
#[non_exhaustive]
pub enum FastqError {
    Io(std::io::Error),
    Format(String),
    Bgzf(String),
    RecordTooLarge { slab_size: usize },
}

impl fmt::Display for FastqError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(err) => write!(f, "I/O error: {err}"),
            Self::Format(msg) => write!(f, "FASTQ parse error: {msg}"),
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
            Self::Format(_) | Self::Bgzf(_) | Self::RecordTooLarge { .. } => None,
        }
    }
}

impl From<std::io::Error> for FastqError {
    fn from(value: std::io::Error) -> Self {
        Self::Io(value)
    }
}
