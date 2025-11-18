use color_eyre::{Result, eyre::Context as _};
use cstr8::CStr8;
use rastair_types::SmallVec;
use rastair_types::SmolStr;
use rust_htslib::bcf::Record;
use std::fmt;

/// A field that can be used in the INFO section.
pub trait FormatField: super::VcfField {
    /// The type of values that this field can hold.
    type Type: FormatFieldValue;

    /// The number of values that can be included with the field.
    const NUMBER: FormatFieldNumber;

    /// The definition of the field for the VCF header.
    fn header() -> String {
        format!(
            "##FORMAT=<ID={},Number={},Type={},Description=\"{}\">",
            Self::ID,
            Self::NUMBER,
            Self::Type::TYPE_NAME,
            Self::DESCRIPTION
        )
    }

    /// Write the field definition to the VCF header.
    fn write_header(header: &mut rust_htslib::bcf::Header) -> Result<()> {
        header.push_record(Self::header().as_bytes());
        Ok(())
    }

    /// Write the field values to the VCF record.
    fn write(&self, record: &mut Record) -> Result<()>;

    /// Description of this field
    fn description() -> Vec<crate::reflect::Format> {
        vec![crate::reflect::Format {
            name: SmolStr::new_static(Self::ID),
            description: SmolStr::new_static(Self::DESCRIPTION),
            number: Self::NUMBER,
            field_type: SmolStr::new_static(Self::Type::TYPE_NAME),
            rust_type: SmolStr::new_static(std::any::type_name::<Self::Type>()),
        }]
    }
}

/// The number of values that can be included with the FORMAT field
///
/// Very simular to [`crate::InfoFieldNumber`], but with additional variants.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FormatFieldNumber {
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
    /// Identical to A except the only alternate alleles defined in the `LAA` field are considered present.
    /// Written as "LA" in the header.
    OnePerLocalAllele,
    /// Identical to R except the only alternate alleles defined in the `LAA` field are considered present.
    /// Written as "LR" in the header.
    OnePerLocalAlleleAndRef,
    /// Identical to G except the only alternate alleles defined in the `LAA` field are considered present.
    /// Written as "LG" in the header.
    OnePerLocalAlleleAndGenotype,
    /// One value for each allele value defined in `GT`.
    /// Written as "P" in the header.
    OnePerAlleleAndGenotype,
    /// One value for each possible base modification for the corresponding ChEBI ID.
    /// Written as "M" in the header.
    OnePerPossibleBaseModification,
}

impl FormatFieldNumber {
    /// Guess the number of values that this field will hold.
    ///
    /// Used for smallvec capacity allocation.
    pub const fn guess_num_values(&self) -> usize {
        match self {
            FormatFieldNumber::Flag => 0,
            FormatFieldNumber::Num(n) => {
                let n = *n as usize;
                if n > 3 { 3 } else { n }
            }
            FormatFieldNumber::OnePerAlt => 1,
            FormatFieldNumber::OnePerAltAndRef => 2,
            FormatFieldNumber::OnePerGenotype => 1,
            FormatFieldNumber::Dot => 1,
            FormatFieldNumber::OnePerLocalAllele => 1,
            FormatFieldNumber::OnePerLocalAlleleAndRef => 2,
            FormatFieldNumber::OnePerLocalAlleleAndGenotype => 1,
            FormatFieldNumber::OnePerAlleleAndGenotype => 1,
            FormatFieldNumber::OnePerPossibleBaseModification => 1,
        }
    }
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Display for FormatFieldNumber {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FormatFieldNumber::Flag => write!(f, "0"),
            FormatFieldNumber::Num(n) => write!(f, "{n}"),
            FormatFieldNumber::OnePerAlt => write!(f, "A"),
            FormatFieldNumber::OnePerAltAndRef => write!(f, "R"),
            FormatFieldNumber::OnePerGenotype => write!(f, "G"),
            FormatFieldNumber::Dot => write!(f, "."),
            FormatFieldNumber::OnePerLocalAllele => write!(f, "LA"),
            FormatFieldNumber::OnePerLocalAlleleAndRef => write!(f, "LR"),
            FormatFieldNumber::OnePerLocalAlleleAndGenotype => write!(f, "LG"),
            FormatFieldNumber::OnePerAlleleAndGenotype => write!(f, "P"),
            FormatFieldNumber::OnePerPossibleBaseModification => write!(f, "M"),
        }
    }
}

/// Types that can be used as values in FORMAT fields.
pub trait FormatFieldValue: Sized {
    /// Possible Types for FORMAT fields are Integer, Float, Character, and String
    const TYPE_NAME: &'static str;

    /// Write the values to the VCF record under the given tag.
    fn write(record: &mut Record, tag: &CStr8, values: &[Self]) -> Result<()>;
}

impl<T: FormatFieldValue + Clone> FormatFieldValue for Option<T> {
    const TYPE_NAME: &'static str = T::TYPE_NAME;

    fn write(record: &mut Record, tag: &CStr8, values: &[Option<T>]) -> Result<()> {
        let non_none_values: Vec<_> = values.iter().filter_map(|v| v.as_ref()).cloned().collect();
        if non_none_values.is_empty() {
            // no values
            Ok(())
        } else {
            T::write(record, tag, non_none_values.as_slice())
        }
    }
}

impl FormatFieldValue for u32 {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &CStr8, values: &[u32]) -> Result<()> {
        record
            .push_format_integer(
                tag,
                &values
                    .iter()
                    .map(|&n| i32::try_from(n))
                    .collect::<Result<SmallVec<i32, 5>, _>>()
                    .wrap_err("Failed to convert u32 to i32")?,
            )
            .wrap_err("Failed to set field")
    }
}

impl FormatFieldValue for u64 {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &CStr8, values: &[u64]) -> Result<()> {
        record
            .push_format_integer(
                tag,
                &values
                    .iter()
                    .map(|&n| i32::try_from(n))
                    .collect::<Result<SmallVec<i32, 5>, _>>()
                    .wrap_err("Failed to convert u64 to i32")?,
            )
            .wrap_err("Failed to set field")
    }
}

impl FormatFieldValue for i32 {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &CStr8, values: &[i32]) -> Result<()> {
        record.push_format_integer(tag, values).wrap_err("Failed to set field")
    }
}

impl FormatFieldValue for i64 {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &CStr8, values: &[i64]) -> Result<()> {
        record
            .push_format_integer(
                tag,
                &values
                    .iter()
                    .map(|&n| i32::try_from(n))
                    .collect::<Result<SmallVec<i32, 5>, _>>()
                    .wrap_err("Failed to convert i64 to i32")?,
            )
            .wrap_err("Failed to set field")
    }
}

impl FormatFieldValue for usize {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &CStr8, values: &[usize]) -> Result<()> {
        record
            .push_format_integer(
                tag,
                &values
                    .iter()
                    .map(|&n| i32::try_from(n))
                    .collect::<Result<SmallVec<i32, 5>, _>>()
                    .wrap_err("Failed to convert usize to i32")?,
            )
            .wrap_err("Failed to set field")
    }
}

impl FormatFieldValue for f32 {
    const TYPE_NAME: &'static str = "Float";

    fn write(record: &mut Record, tag: &CStr8, values: &[f32]) -> Result<()> {
        record.push_format_float(tag, values).wrap_err("Failed to set field")
    }
}

impl FormatFieldValue for f64 {
    const TYPE_NAME: &'static str = "Float";

    #[allow(clippy::cast_possible_truncation)] // Allow casting f64 to f32, which is common in VCF
    fn write(record: &mut Record, tag: &CStr8, values: &[f64]) -> Result<()> {
        record
            .push_format_float(tag, &values.iter().map(|&n| n as f32).collect::<SmallVec<f32, 5>>())
            .wrap_err("Failed to set field")
    }
}

impl FormatFieldValue for String {
    const TYPE_NAME: &'static str = "String";

    fn write(record: &mut Record, tag: &CStr8, values: &[String]) -> Result<()> {
        record
            .push_format_string(
                tag,
                &values.iter().map(|s| s.as_bytes()).collect::<SmallVec<&[u8], 5>>(),
            )
            .wrap_err("Failed to set field")
    }
}

impl FormatFieldValue for SmolStr {
    const TYPE_NAME: &'static str = "String";

    fn write(record: &mut Record, tag: &CStr8, values: &[SmolStr]) -> Result<()> {
        record
            .push_format_string(
                tag,
                &values.iter().map(|s| s.as_bytes()).collect::<SmallVec<&[u8], 5>>(),
            )
            .wrap_err("Failed to set field")
    }
}

/// Define a VCF format field.
///
/// # Syntax
///
/// ```rust
/// use rastair_vcf::{format_field, FormatFieldNumber};
/// type Type = u32; // or any other type that implements FormatFieldValue
///
/// format_field!(Name(Type), "ID", "Description", FormatFieldNumber::OnePerAlt);
/// ```
///
/// This will define a struct `Name` that implements the [`FormatField`] trait
/// (as well as [`crate::VcfField`] and [`crate::HeaderField`]).
///
/// The last parameter must be a variant of [`FormatFieldNumber`].
///
/// # Inlining
///
/// Fields can be lists and we want to keep most of the data on the
/// stack/inline. The generated struct will use [`SmallVec`](smallvec::SmallVec)
/// to inline a guessed number of values baed on the given [`FormatFieldNumber`].
#[macro_export]
macro_rules! format_field {
    ($name:ident($type:ty), $id:expr, $desc:expr, 1) => {
        #[doc = $desc]
        #[doc = "("]
        #[doc = stringify!($number)]
        #[doc = ")"]
        #[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
        #[repr(transparent)]
        pub struct $name(pub $type);

        impl std::ops::Deref for $name {
            type Target = $type;

            fn deref(&self) -> &Self::Target {
                &self.0
            }
        }

        impl Default for $name {
            fn default() -> Self {
                Self(Default::default())
            }
        }

        impl $crate::VcfField for $name {
            const ID: &'static cstr8::CStr8 = cstr8::cstr8!($id);
        }

        impl $crate::HeaderField for $name {
            const DESCRIPTION: &'static str = $desc;
        }

        impl $crate::FormatField for $name {
            type Type = $type;
            const NUMBER: $crate::FormatFieldNumber = $crate::FormatFieldNumber::Num(1);

            fn write(&self, record: &mut rust_htslib::bcf::Record) -> color_eyre::Result<()> {
                use color_eyre::eyre::WrapErr;
                use $crate::{FormatFieldValue, VcfField as _};

                <$type as $crate::FormatFieldValue>::write(record, Self::ID, &[self.0])
                    .wrap_err_with(|| {
                        format!(
                            "Failed to write format field {} (type {})",
                            Self::ID,
                            <Self::Type as FormatFieldValue>::TYPE_NAME
                        )
                    })
            }
        }
    };

    ($name:ident($type:ty), $id:expr, $desc:expr, $number:expr) => {
        #[doc = $desc]
        #[doc = "("]
        #[doc = stringify!($number)]
        #[doc = ")"]
        #[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
        #[repr(transparent)]
        pub struct $name(pub rastair_types::SmallVec<$type, { $number.guess_num_values() }>);

        impl std::ops::Deref for $name {
            type Target = [$type];

            fn deref(&self) -> &Self::Target {
                &self.0.as_slice()
            }
        }

        impl Default for $name {
            fn default() -> Self {
                Self(Default::default())
            }
        }

        impl $crate::VcfField for $name {
            const ID: &'static cstr8::CStr8 = cstr8::cstr8!($id);
        }

        impl $crate::HeaderField for $name {
            const DESCRIPTION: &'static str = $desc;
        }

        impl $crate::FormatField for $name {
            type Type = $type;
            const NUMBER: $crate::FormatFieldNumber = $number;

            fn write(&self, record: &mut rust_htslib::bcf::Record) -> color_eyre::Result<()> {
                use color_eyre::eyre::WrapErr;
                use $crate::{FormatFieldValue, VcfField as _};

                <$type as $crate::FormatFieldValue>::write(record, Self::ID, &self.0).wrap_err_with(
                    || {
                        format!(
                            "Failed to write format field {} (type {})",
                            Self::ID,
                            <Self::Type as FormatFieldValue>::TYPE_NAME
                        )
                    },
                )
            }
        }
    };
}

#[cfg(test)]
mod tests {
    use super::super::FormatFieldNumber;
    use super::*;
    use insta::assert_snapshot;
    use rastair_types::smallvec::smallvec;
    use rust_htslib::bcf::{Format, Header, Writer};
    use tempfile::TempDir;

    #[test]
    fn format_header() {
        format_field!(Foo(String), "GT", "Foo", FormatFieldNumber::OnePerGenotype);

        assert_snapshot!(
            Foo::header(),
            @r###"##FORMAT=<ID=GT,Number=G,Type=String,Description="Foo">"###
        );
    }

    #[test]
    fn integers() -> Result<()> {
        format_field!(FieldU32(u32), "U32", "Test u32", FormatFieldNumber::OnePerAlt);
        format_field!(FieldU64(u64), "U64", "Test u64", FormatFieldNumber::OnePerAlt);
        format_field!(FieldI32(i32), "I32", "Test i32", FormatFieldNumber::OnePerAlt);
        format_field!(FieldI64(i64), "I64", "Test i64", FormatFieldNumber::OnePerAlt);
        format_field!(FieldUsize(usize), "Usize", "Test usize", FormatFieldNumber::OnePerAlt);

        let mut header = Header::new();
        header.push_record(b"##fileformat=VCFv4.2");
        header.push_record(br#"##contig=<ID=1,length=10>"#);
        header.push_record(FieldU32::header().as_bytes());
        header.push_record(FieldU64::header().as_bytes());
        header.push_record(FieldI32::header().as_bytes());
        header.push_record(FieldI64::header().as_bytes());
        header.push_record(FieldUsize::header().as_bytes());
        header.push_sample(b"one");

        let temp_dir = TempDir::new()?;
        let temp_file = temp_dir.path().join("test.vcf");
        let mut vcf = Writer::from_path(&temp_file, &header, true, Format::Vcf)?;
        let mut record = vcf.empty_record();

        // different types
        FieldU32(smallvec![1]).write(&mut record)?;
        assert_eq!((*record.format(b"U32").integer()?)[0], &[1i32]);

        FieldI32(smallvec![42]).write(&mut record)?;
        assert_eq!((*record.format(b"I32").integer()?)[0], &[42i32]);

        FieldU64(smallvec![42]).write(&mut record)?;
        assert_eq!((*record.format(b"U64").integer()?)[0], &[42i32]);

        FieldI64(smallvec![42]).write(&mut record)?;
        assert_eq!((*record.format(b"I64").integer()?)[0], &[42i32]);

        FieldUsize(smallvec![42]).write(&mut record)?;
        assert_eq!((*record.format(b"Usize").integer()?)[0], &[42i32]);

        // lists
        FieldU32(smallvec![1, 2]).write(&mut record)?;
        assert_eq!((*record.format(b"U32").integer()?)[0], &[1i32, 2]);

        FieldI64(smallvec![1, 2]).write(&mut record)?;
        assert_eq!((*record.format(b"I64").integer()?)[0], &[1i32, 2]);

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
}
