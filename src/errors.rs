use thiserror::Error;

#[derive(Debug, Error)]
pub enum HorScanError {
    #[error("I/O error: {0}")]
    Io(String),
    #[error("Parse error: {0}")]
    ParseError(String),
    #[error("Invalid data: {0}")]
    InvalidData(String),
}
