use thiserror::Error;

#[derive(Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Base {
    A = b'A',
    C = b'C',
    G = b'G',
    T = b'T',
}

impl std::fmt::Display for Base {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", (*self) as u8 as char)
    }
}

impl Base {
    pub fn display_colored(&self) -> &str {
        match self {
            Base::A => "\x1b[32mA\x1b[0m", // green
            Base::C => "\x1b[34mC\x1b[0m", // blue
            Base::G => "\x1b[33mG\x1b[0m", // yellow
            Base::T => "\x1b[31mT\x1b[0m", // red
        }
    }
}

impl std::ops::Deref for Base {
    type Target = u8;

    fn deref(&self) -> &Self::Target {
        match self {
            Base::A => &b'A',
            Base::C => &b'C',
            Base::G => &b'G',
            Base::T => &b'T',
        }
    }
}

impl std::fmt::Debug for Base {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self)
    }
}

impl std::str::FromStr for Base {
    type Err = BaseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();
        let Some(first) = s.as_bytes().get(0) else {
            return Err(BaseError::Empty);
        };
        first.as_base()
    }
}

#[derive(Debug, Error)]
pub enum BaseError {
    #[error("Empty")]
    Empty,
    #[error("Invalid base `{0}`")]
    InvalidBaseError(u8),
}

pub trait TryAsBase {
    fn as_base(&self) -> Result<Base, BaseError>;
}

impl TryAsBase for u8 {
    fn as_base(&self) -> Result<Base, BaseError> {
        match self {
            b'A' => Ok(Base::A),
            b'C' => Ok(Base::C),
            b'G' => Ok(Base::G),
            b'T' => Ok(Base::T),
            _ => Err(BaseError::InvalidBaseError(*self)),
        }
    }
}
