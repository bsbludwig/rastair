use smol_str::SmolStr;
use thiserror::Error;

#[derive(Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[repr(u8)]
pub enum Base {
    A = b'A',
    C = b'C',
    G = b'G',
    T = b'T',
    Unknown = b'N',
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
            Base::Unknown => "N",          // white
        }
    }

    /// Get the inverse base (complementary base)
    pub fn inverse(&self) -> Base {
        match self {
            Base::A => Base::T,
            Base::C => Base::G,
            Base::G => Base::C,
            Base::T => Base::A,
            Base::Unknown => Base::Unknown,
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Base::A => "A",
            Base::C => "C",
            Base::G => "G",
            Base::T => "T",
            Base::Unknown => "N",
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
            Base::Unknown => &b'N',
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
        Ok(Base::from(*first))
    }
}

impl From<u8> for Base {
    fn from(value: u8) -> Self {
        match value {
            b'a' => Base::A,
            b'c' => Base::C,
            b'g' => Base::G,
            b't' => Base::T,
            b'A' => Base::A,
            b'C' => Base::C,
            b'G' => Base::G,
            b'T' => Base::T,
            _ => Base::Unknown,
        }
    }
}

impl From<&u8> for Base {
    fn from(value: &u8) -> Self {
        Base::from(*value)
    }
}

impl From<Base> for SmolStr {
    fn from(val: Base) -> Self {
        SmolStr::new_inline(val.into())
    }
}

impl From<Base> for &'static str {
    fn from(val: Base) -> Self {
        val.as_str()
    }
}

#[derive(Debug, Error)]
pub enum BaseError {
    #[error("Empty")]
    Empty,
    #[error("Invalid base {base}", base=(0 as char))]
    InvalidBaseError(u8),
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    proptest::proptest! {
        #[test]
        fn proptest_roundtrip(input: u8) {
            let base = match Base::from(input) {
                Base::Unknown => return Ok(()), // We don't panic on unknown bases
                base => base
            };
            assert_eq!(*base, (input as char).to_ascii_uppercase() as u8);
        }

        #[test]
        fn proptest_roundtrip_str(input in r"\PC{0,10}" ) {
            let base = match Base::from_str(&input) {
                Err(_) => return Ok(()), // We don't panic on invalid bases
                Ok(Base::Unknown) => return Ok(()), // We don't panic on unknown bases
                Ok(base) => base,
            };
            assert_eq!(*base, input.trim().to_ascii_uppercase().as_bytes()[0]);
        }
    }

    #[test]
    fn test_u8_to_base_valid() {
        let valid_bases = [b'A', b'C', b'G', b'T'];
        for &base in &valid_bases {
            let parsed = Base::from(base);
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
