use smol_str::SmolStr;
use thiserror::Error;

#[derive(Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Base {
    A = b'A',
    C = b'C',
    G = b'G',
    T = b'T',
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl std::fmt::Display for Base {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", (*self) as u8 as char)
    }
}

impl Base {
    #[cfg_attr(coverage_nightly, coverage(off))]
    pub fn display_colored(&self) -> &str {
        match self {
            Base::A => "\x1b[32mA\x1b[0m", // green
            Base::C => "\x1b[34mC\x1b[0m", // blue
            Base::G => "\x1b[33mG\x1b[0m", // yellow
            Base::T => "\x1b[31mT\x1b[0m", // red
        }
    }

    /// Get the inverse base (complementary base)
    pub fn inverse(&self) -> Base {
        match self {
            Base::A => Base::T,
            Base::C => Base::G,
            Base::G => Base::C,
            Base::T => Base::A,
        }
    }

    pub fn as_char(&self) -> char {
        (*self) as u8 as char
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

#[cfg_attr(coverage_nightly, coverage(off))]
impl std::fmt::Debug for Base {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if f.alternate() || cfg!(test) {
            write!(f, "{}", (*self) as u8 as char)
        } else {
            write!(f, "{}", self.display_colored())
        }
    }
}

impl std::str::FromStr for Base {
    type Err = BaseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();
        let Some(first) = s.as_bytes().first() else {
            return Err(BaseError::Empty);
        };
        first.as_base()
    }
}

impl TryFrom<u8> for Base {
    type Error = BaseError;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            b'a' => Ok(Base::A),
            b'c' => Ok(Base::C),
            b'g' => Ok(Base::G),
            b't' => Ok(Base::T),
            b'A' => Ok(Base::A),
            b'C' => Ok(Base::C),
            b'G' => Ok(Base::G),
            b'T' => Ok(Base::T),
            _ => Err(BaseError::InvalidBaseError(value)),
        }
    }
}

impl From<Base> for SmolStr {
    fn from(val: Base) -> Self {
        SmolStr::new_inline(val.into())
    }
}

impl From<Base> for &'static str {
    fn from(val: Base) -> Self {
        match val {
            Base::A => "A",
            Base::C => "C",
            Base::G => "G",
            Base::T => "T",
        }
    }
}

#[derive(Debug, Error)]
pub enum BaseError {
    #[error("Empty")]
    Empty,
    #[error("Invalid base {base}", base=(0 as char))]
    InvalidBaseError(u8),
}

pub trait TryAsBase {
    fn as_base(&self) -> Result<Base, BaseError>;
}

impl TryAsBase for u8 {
    fn as_base(&self) -> Result<Base, BaseError> {
        (*self).try_into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    proptest::proptest! {
        #[test]
        fn proptest_roundtrip(input: u8) {
            let Ok(base) = input.as_base() else {
                // We're just checking that there is no panic, but errors are fine!
                return Ok(());
            };
            assert_eq!(*base, (input as char).to_ascii_uppercase() as u8);
        }

        #[test]
        fn proptest_roundtrip_str(input in r"\PC{0,10}" ) {
            let Ok(base) = Base::from_str(&input) else {
                // We're just checking that there is no panic, but errors are fine!
                return Ok(());
            };
            assert_eq!(*base, input.trim().to_ascii_uppercase().as_bytes()[0]);
        }
    }

    #[test]
    fn test_u8_to_base_valid() {
        let valid_bases = [b'A', b'C', b'G', b'T'];
        for &base in &valid_bases {
            let parsed = base.as_base().unwrap();
            assert_eq!(*parsed, base);
            // display
            assert_eq!(parsed.to_string(), (base as char).to_string());
            // debug -- same, actually
            assert_eq!(format!("{parsed:#?}"), (base as char).to_string());
            // display in color
            let colored = parsed.display_colored();
            // skip initial ANSI color codes
            let x = colored.chars().nth(5).unwrap();
            assert_eq!(x, base as char);
        }
    }
}
