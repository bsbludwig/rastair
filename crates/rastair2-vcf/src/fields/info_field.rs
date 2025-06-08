use std::fmt;

use color_eyre::{Result, eyre::Context as _};
use rust_htslib::bcf::Record;
use smallvec::SmallVec;
use smol_str::SmolStr;

/// A field that can be used in the INFO section.
pub trait InfoField: super::VcfField {
    /// The type of values that this field can hold.
    type Type: InfoFieldValue;

    /// The number of values that can be included with the field.
    const NUMBER: InfoFieldNumber;

    /// The definition of the field for the VCF header.
    fn header() -> String {
        format!(
            "##INFO=<ID={},Number={},Type={},Description=\"{}\">",
            Self::ID,
            Self::NUMBER,
            Self::Type::TYPE_NAME,
            Self::DESCRIPTION
        )
    }

    /// Write the field values to the VCF record.
    fn write(&self, record: &mut Record) -> Result<()>;
}

/// The number of values that can be included with the INFO field
pub enum InfoFieldNumber {
    /// A flag, no values. Written as "0" in the header.
    Flag,
    /// A fixed number of values, must be non-zero. Written as a number in the header.
    Num(u32),
    /// One value per alternative allele. Written as "A" in the header.
    OnePerAlt,
    /// One value per alternative allele and reference allele. Written as "R" in the header.
    OnePerAltAndRef,
    /// One value per genotype. Written as "G" in the header.
    OnePerGenotype,
    /// Variable number of values or unknown or unbounded. Written as "." in the header.
    Dot,
}

impl InfoFieldNumber {
    /// Guess the number of values that this field will hold.
    ///
    /// Used for smallvec capacity allocation.
    pub const fn guess_num_values(&self) -> usize {
        match self {
            InfoFieldNumber::Flag => 0,
            InfoFieldNumber::Num(n) => {
                let n = *n as usize;
                if n > 3 { 3 } else { n }
            }
            InfoFieldNumber::OnePerAlt => 1,
            InfoFieldNumber::OnePerAltAndRef => 3,
            InfoFieldNumber::OnePerGenotype => 2, // This is per genotype, not per sample
            InfoFieldNumber::Dot => 1,            // Represents variable or unknown number of values
        }
    }
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Display for InfoFieldNumber {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InfoFieldNumber::Flag => write!(f, "0"),
            InfoFieldNumber::Num(n) => write!(f, "{n}"),
            InfoFieldNumber::OnePerAlt => write!(f, "A"),
            InfoFieldNumber::OnePerAltAndRef => write!(f, "R"),
            InfoFieldNumber::OnePerGenotype => write!(f, "G"),
            InfoFieldNumber::Dot => write!(f, "."),
        }
    }
}

/// Types that can be used as values in INFO fields.
pub trait InfoFieldValue: Sized {
    /// Possible Types for INFO fields are: Integer, Float, Flag, Character, and String
    const TYPE_NAME: &'static str;

    /// Write the values to the VCF record under the given tag.
    fn write(record: &mut Record, tag: &str, values: &[Self]) -> Result<()>;
}

impl InfoFieldValue for () {
    const TYPE_NAME: &'static str = "Flag";

    fn write(record: &mut Record, tag: &str, _values: &[()]) -> Result<()> {
        record
            .clear_info_flag(tag.as_bytes())
            .wrap_err_with(|| format!("Failed to clear info field {tag} (Flag)"))?;
        record.push_info_flag(tag.as_bytes()).wrap_err_with(|| format!("Failed to set flag {tag}"))
    }
}

impl InfoFieldValue for u32 {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &str, values: &[u32]) -> Result<()> {
        record
            .clear_info_integer(tag.as_bytes())
            .wrap_err_with(|| format!("Failed to clear info field {tag} ({})", Self::TYPE_NAME))?;
        record
            .push_info_integer(
                tag.as_bytes(),
                &values
                    .iter()
                    .map(|&n| i32::try_from(n))
                    .collect::<Result<SmallVec<i32, 5>, _>>()
                    .wrap_err_with(|| {
                        format!("Failed to convert u32 to i32 for info field {tag}")
                    })?,
            )
            .wrap_err_with(|| format!("Failed to set info field {tag} ({})", Self::TYPE_NAME))
    }
}

impl InfoFieldValue for u64 {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &str, values: &[u64]) -> Result<()> {
        record
            .clear_info_integer(tag.as_bytes())
            .wrap_err_with(|| format!("Failed to clear info field {tag} ({})", Self::TYPE_NAME))?;
        record
            .push_info_integer(
                tag.as_bytes(),
                &values
                    .iter()
                    .map(|&n| i32::try_from(n))
                    .collect::<Result<SmallVec<i32, 5>, _>>()
                    .wrap_err_with(|| {
                        format!("Failed to convert u64 to i32 for info field {tag}")
                    })?,
            )
            .wrap_err_with(|| format!("Failed to set info field {tag} ({})", Self::TYPE_NAME))
    }
}

impl InfoFieldValue for i32 {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &str, values: &[i32]) -> Result<()> {
        record
            .clear_info_integer(tag.as_bytes())
            .wrap_err_with(|| format!("Failed to clear info field {tag} ({})", Self::TYPE_NAME))?;
        record
            .push_info_integer(tag.as_bytes(), values)
            .wrap_err_with(|| format!("Failed to set info field {tag} ({})", Self::TYPE_NAME))
    }
}

impl InfoFieldValue for i64 {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &str, values: &[i64]) -> Result<()> {
        record
            .clear_info_integer(tag.as_bytes())
            .wrap_err_with(|| format!("Failed to clear info field {tag} ({})", Self::TYPE_NAME))?;
        record
            .push_info_integer(
                tag.as_bytes(),
                &values
                    .iter()
                    .map(|&n| i32::try_from(n))
                    .collect::<Result<SmallVec<i32, 5>, _>>()
                    .wrap_err_with(|| {
                        format!("Failed to convert i64 to i32 for info field {tag}")
                    })?,
            )
            .wrap_err_with(|| format!("Failed to set info field {tag} ({})", Self::TYPE_NAME))
    }
}

impl InfoFieldValue for usize {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &str, values: &[usize]) -> Result<()> {
        record
            .clear_info_integer(tag.as_bytes())
            .wrap_err_with(|| format!("Failed to clear info field {tag} ({})", Self::TYPE_NAME))?;
        record
            .push_info_integer(
                tag.as_bytes(),
                &values
                    .iter()
                    .map(|&n| i32::try_from(n))
                    .collect::<Result<SmallVec<i32, 5>, _>>()
                    .wrap_err_with(|| {
                        format!("Failed to convert usize to i32 for info field {tag}")
                    })?,
            )
            .wrap_err_with(|| format!("Failed to set info field {tag} ({})", Self::TYPE_NAME))
    }
}

impl InfoFieldValue for f32 {
    const TYPE_NAME: &'static str = "Float";

    fn write(record: &mut Record, tag: &str, values: &[f32]) -> Result<()> {
        record
            .clear_info_float(tag.as_bytes())
            .wrap_err_with(|| format!("Failed to clear info field {tag} ({})", Self::TYPE_NAME))?;
        record
            .push_info_float(tag.as_bytes(), values)
            .wrap_err_with(|| format!("Failed to set info field {tag} ({})", Self::TYPE_NAME))
    }
}

impl InfoFieldValue for f64 {
    const TYPE_NAME: &'static str = "Float";

    #[allow(clippy::cast_possible_truncation)] // Allow casting f64 to f32, which is common in VCF
    fn write(record: &mut Record, tag: &str, values: &[f64]) -> Result<()> {
        record
            .clear_info_float(tag.as_bytes())
            .wrap_err_with(|| format!("Failed to clear info field {tag} ({})", Self::TYPE_NAME))?;
        record
            .push_info_float(
                tag.as_bytes(),
                &values.iter().map(|&n| n as f32).collect::<SmallVec<f32, 5>>(),
            )
            .wrap_err_with(|| format!("Failed to set info field {tag} ({})", Self::TYPE_NAME))
    }
}

impl InfoFieldValue for String {
    const TYPE_NAME: &'static str = "String";

    fn write(record: &mut Record, tag: &str, values: &[String]) -> Result<()> {
        record
            .push_info_string(
                tag.as_bytes(),
                &values.iter().map(|s| s.as_bytes()).collect::<SmallVec<&[u8], 5>>(),
            )
            .wrap_err_with(|| format!("Failed to set info field {tag} ({})", Self::TYPE_NAME))
    }
}

impl InfoFieldValue for SmolStr {
    const TYPE_NAME: &'static str = "String";

    fn write(record: &mut Record, tag: &str, values: &[SmolStr]) -> Result<()> {
        record
            .push_info_string(
                tag.as_bytes(),
                &values.iter().map(|s| s.as_bytes()).collect::<SmallVec<&[u8], 5>>(),
            )
            .wrap_err_with(|| format!("Failed to set info field {tag} ({})", Self::TYPE_NAME))
    }
}

/// Define a VCF info field.
///
/// # Syntax
///
/// ```rust
/// use rastair2_vcf::{info_field, InfoFieldNumber};
/// # type Type = u32; // or any other type that implements InfoFieldValue
///
/// info_field!(Name(Type), "ID", "Description", InfoFieldNumber::OnePerAlt);
/// ```
///
/// This will define a struct `Name` that implements the [`InfoField`] trait (as
/// well as [`crate::VcfField`] and [`crate::HeaderField`]).
///
/// The last parameter must be a variant of [`InfoFieldNumber`].
///
/// # Inlining
///
/// Fields can be lists and we want to keep most of the data on the
/// stack/inline. The generated struct will use [`SmallVec`](smallvec::SmallVec)
/// to inline a guessed number of values baed on the given [`InfoFieldNumber`].
#[macro_export]
macro_rules! info_field {
    ($name:ident($type:tt), $id:expr, $desc:expr, $number:expr) => {
        #[doc = $desc]
        #[doc = "("]
        #[doc = stringify!($number)]
        #[doc = ")"]
        #[derive(Debug, Clone)]
        pub struct $name(pub smallvec::SmallVec<$type, { $number.guess_num_values() }>);

        impl std::ops::Deref for $name {
            type Target = [$type];

            fn deref(&self) -> &Self::Target {
                &self.0.as_slice()
            }
        }

        impl $crate::VcfField for $name {
            const ID: &'static str = $id;
        }

        impl $crate::HeaderField for $name {
            const DESCRIPTION: &'static str = $desc;
        }

        impl $crate::InfoField for $name {
            type Type = $type;
            const NUMBER: $crate::InfoFieldNumber = $number;

            fn write(&self, record: &mut rust_htslib::bcf::Record) -> color_eyre::Result<()> {
                use $crate::VcfField as _;

                <$type as $crate::InfoFieldValue>::write(record, Self::ID, &self.0)
            }
        }
    };

    ($name:ident, $id:expr, $desc:expr) => {
        #[doc = $desc]
        #[doc = "(flag)"]
        #[derive(Debug, Clone)]
        pub struct $name;

        impl $crate::VcfField for $name {
            const ID: &'static str = $id;
        }

        impl $crate::HeaderField for $name {
            const DESCRIPTION: &'static str = $desc;
        }

        impl $crate::InfoField for $name {
            type Type = ();
            const NUMBER: $crate::InfoFieldNumber = $crate::InfoFieldNumber::Flag;

            fn write(&self, record: &mut rust_htslib::bcf::Record) -> color_eyre::Result<()> {
                use $crate::VcfField as _;

                <() as $crate::InfoFieldValue>::write(record, Self::ID, &[])
            }
        }
    };
}

#[cfg(test)]
mod tests {
    use super::super::{InfoFieldNumber, VcfField};
    use super::*;
    use color_eyre::{Result, eyre::ContextCompat};
    use insta::assert_snapshot;
    use rust_htslib::bcf::{Format, Header, Writer};
    use smallvec::smallvec;
    use tempfile::TempDir;

    #[test]
    fn info_header() {
        info_field!(AlleleFrequency(f64), "AF", "Allele Frequency", InfoFieldNumber::OnePerAlt);

        assert_snapshot!(
            AlleleFrequency::header(),
            @r###"##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">"###
        );
    }

    #[test]
    fn flags() -> Result<()> {
        info_field!(Flag, "Flag", "Test flag 1");
        info_field!(Glag, "Glag", "Test flag 2");
        info_field!(Klag, "Klag", "Test flag 3");

        let mut header = Header::new();
        header.push_record(b"##fileformat=VCFv4.2");
        header.push_record(br#"##contig=<ID=1,length=10>"#);
        header.push_record(Flag::header().as_bytes());
        header.push_record(Glag::header().as_bytes());
        header.push_record(Klag::header().as_bytes());

        let temp_dir = TempDir::new()?;
        let temp_file = temp_dir.path().join("test.vcf");
        let mut vcf = Writer::from_path(&temp_file, &header, true, Format::Vcf)?;
        let mut record = vcf.empty_record();

        assert!(!record.info(Flag::ID.as_bytes()).flag()?);
        Flag.write(&mut record)?;
        assert!(record.info(Flag::ID.as_bytes()).flag()?);

        Glag.write(&mut record)?;
        assert!(record.info(b"Glag").flag()?);

        Klag.write(&mut record)?;
        assert!(record.info(b"Klag").flag()?);

        vcf.write(&record).wrap_err("Failed to write record")?;
        drop(vcf);
        let result = std::fs::read_to_string(temp_file).wrap_err("read temp file")?;
        assert_snapshot!(result);

        Ok(())
    }

    #[test]
    fn integers() -> Result<()> {
        info_field!(FieldU32(u32), "U32", "Test u32", InfoFieldNumber::OnePerAlt);
        info_field!(FieldU64(u64), "U64", "Test u64", InfoFieldNumber::OnePerAlt);
        info_field!(FieldI32(i32), "I32", "Test i32", InfoFieldNumber::OnePerAlt);
        info_field!(FieldI64(i64), "I64", "Test i64", InfoFieldNumber::OnePerAlt);
        info_field!(FieldUsize(usize), "Usize", "Test usize", InfoFieldNumber::OnePerAlt);

        let mut header = Header::new();
        header.push_record(b"##fileformat=VCFv4.2");
        header.push_record(br#"##contig=<ID=1,length=10>"#);
        header.push_record(FieldU32::header().as_bytes());
        header.push_record(FieldU64::header().as_bytes());
        header.push_record(FieldI32::header().as_bytes());
        header.push_record(FieldI64::header().as_bytes());
        header.push_record(FieldUsize::header().as_bytes());

        let temp_dir = TempDir::new()?;
        let temp_file = temp_dir.path().join("test.vcf");
        let mut vcf = Writer::from_path(&temp_file, &header, true, Format::Vcf)?;
        let mut record = vcf.empty_record();

        // different types
        FieldU32(smallvec![1]).write(&mut record)?;
        assert_eq!(*record.info(b"U32").integer()?.wrap_err("none")?, &[1i32]);

        FieldI32(smallvec![42]).write(&mut record)?;
        assert_eq!(*record.info(b"I32").integer()?.wrap_err("none")?, &[42i32]);

        FieldU64(smallvec![42]).write(&mut record)?;
        assert_eq!(*record.info(b"U64").integer()?.wrap_err("none")?, &[42i32]);

        FieldI64(smallvec![42]).write(&mut record)?;
        assert_eq!(*record.info(b"I64").integer()?.wrap_err("none")?, &[42i32]);

        FieldUsize(smallvec![42]).write(&mut record)?;
        assert_eq!(*record.info(b"Usize").integer()?.wrap_err("none")?, &[42i32]);

        // lists
        FieldU32(smallvec![1, 2]).write(&mut record)?;
        assert_eq!(*record.info(b"U32").integer()?.wrap_err("none")?, &[1i32, 2]);

        FieldI64(smallvec![1, 2]).write(&mut record)?;
        assert_eq!(*record.info(b"I64").integer()?.wrap_err("none")?, &[1i32, 2]);

        vcf.write(&record).wrap_err("Failed to write record")?;
        drop(vcf);
        let result = std::fs::read_to_string(temp_file).wrap_err("read temp file")?;
        assert_snapshot!(result);

        // test that we catch overflow
        assert!(FieldU32(smallvec![u32::MAX]).write(&mut record).is_err());
        assert!(FieldU64(smallvec![u64::MAX]).write(&mut record).is_err());
        assert!(FieldI64(smallvec![i64::MAX]).write(&mut record).is_err());
        // but i32 is base type so it's fine
        assert!(FieldI32(smallvec![i32::MAX]).write(&mut record).is_ok());

        Ok(())
    }

    #[test]
    fn floats() -> Result<()> {
        info_field!(FieldF32(f32), "F32", "Test f32", InfoFieldNumber::OnePerAlt);
        info_field!(FieldF64(f64), "F64", "Test f64", InfoFieldNumber::OnePerAlt);

        let mut header = Header::new();
        header.push_record(b"##fileformat=VCFv4.2");
        header.push_record(br#"##contig=<ID=1,length=10>"#);
        header.push_record(FieldF32::header().as_bytes());
        header.push_record(FieldF64::header().as_bytes());

        let temp_dir = TempDir::new()?;
        let temp_file = temp_dir.path().join("test.vcf");
        let mut vcf = Writer::from_path(&temp_file, &header, true, Format::Vcf)?;
        let mut record = vcf.empty_record();

        FieldF32(smallvec![1.4]).write(&mut record)?;
        assert_eq!(*record.info(b"F32").float()?.wrap_err("none")?, &[1.4f32]);
        FieldF32(smallvec![1.1, 2.2]).write(&mut record)?;
        assert_eq!(*record.info(b"F32").float()?.wrap_err("none")?, &[1.1f32, 2.2]);

        FieldF64(smallvec![42.42]).write(&mut record)?;
        assert_eq!(*record.info(b"F64").float()?.wrap_err("none")?, &[42.42f32]);

        vcf.write(&record).wrap_err("Failed to write record")?;
        drop(vcf);
        let result = std::fs::read_to_string(temp_file).wrap_err("read temp file")?;
        assert_snapshot!(result);

        Ok(())
    }

    #[test]
    fn strings() -> Result<()> {
        info_field!(FieldString(String), "STR", "Test string", InfoFieldNumber::OnePerAlt);

        let mut header = Header::new();
        header.push_record(b"##fileformat=VCFv4.2");
        header.push_record(br#"##contig=<ID=1,length=10>"#);
        header.push_record(FieldString::header().as_bytes());

        let temp_dir = TempDir::new()?;
        let temp_file = temp_dir.path().join("test.vcf");
        let mut vcf = Writer::from_path(&temp_file, &header, true, Format::Vcf)?;
        let mut record = vcf.empty_record();

        FieldString(smallvec!["test".into()]).write(&mut record)?;
        assert_eq!(*record.info(b"STR").string()?.wrap_err("none")?, &[b"test"]);

        FieldString(smallvec!["test1".into(), "test2".into()]).write(&mut record)?;
        assert_eq!(*record.info(b"STR").string()?.wrap_err("none")?, &[b"test1", b"test2"]);

        vcf.write(&record).wrap_err("Failed to write record")?;
        drop(vcf);
        let result = std::fs::read_to_string(temp_file).wrap_err("read temp file")?;
        assert_snapshot!(result);

        Ok(())
    }
}
