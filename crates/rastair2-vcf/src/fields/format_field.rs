use color_eyre::{Result, eyre::Context as _};
use rust_htslib::bcf::Record;
use smallvec::SmallVec;
use std::fmt::Display;

/// A field that can be used in the INFO section.
pub trait FormatField: super::VcfField {
    /// The type of values that this field can hold.
    type Type: FormatFieldValue + Display;

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

    /// Write the field values to the VCF record.
    fn write(&self, record: &mut Record) -> Result<()>;
}

/// Types that can be used as values in FORMAT fields.
pub trait FormatFieldValue: Sized {
    /// Possible Types for FORMAT fields are
    // TODO: fill in
    const TYPE_NAME: &'static str;

    /// Write the values to the VCF record under the given tag.
    fn write(record: &mut Record, tag: &str, values: &[Self]) -> Result<()>;
}

impl FormatFieldValue for u32 {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &str, values: &[u32]) -> Result<()> {
        record
            .push_format_integer(
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

impl FormatFieldValue for u64 {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &str, values: &[u64]) -> Result<()> {
        record
            .push_format_integer(
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

impl FormatFieldValue for i32 {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &str, values: &[i32]) -> Result<()> {
        record
            .push_format_integer(tag.as_bytes(), values)
            .wrap_err_with(|| format!("Failed to set info field {tag} ({})", Self::TYPE_NAME))
    }
}

impl FormatFieldValue for i64 {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &str, values: &[i64]) -> Result<()> {
        record
            .push_format_integer(
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

impl FormatFieldValue for usize {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &str, values: &[usize]) -> Result<()> {
        record
            .push_format_integer(
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

impl FormatFieldValue for f32 {
    const TYPE_NAME: &'static str = "Float";

    fn write(record: &mut Record, tag: &str, values: &[f32]) -> Result<()> {
        record
            .push_format_float(tag.as_bytes(), values)
            .wrap_err_with(|| format!("Failed to set info field {tag} ({})", Self::TYPE_NAME))
    }
}

impl FormatFieldValue for f64 {
    const TYPE_NAME: &'static str = "Float";

    #[allow(clippy::cast_possible_truncation)] // Allow casting f64 to f32, which is common in VCF
    fn write(record: &mut Record, tag: &str, values: &[f64]) -> Result<()> {
        record
            .push_format_float(
                tag.as_bytes(),
                &values.iter().map(|&n| n as f32).collect::<SmallVec<f32, 5>>(),
            )
            .wrap_err_with(|| format!("Failed to set info field {tag} ({})", Self::TYPE_NAME))
    }
}

impl FormatFieldValue for String {
    const TYPE_NAME: &'static str = "String";

    fn write(record: &mut Record, tag: &str, values: &[String]) -> Result<()> {
        record
            .push_format_string(
                tag.as_bytes(),
                &values.iter().map(|s| s.as_bytes()).collect::<SmallVec<&[u8], 5>>(),
            )
            .wrap_err_with(|| format!("Failed to set info field {tag} ({})", Self::TYPE_NAME))
    }
}

/// Define a VCF format field.
///
/// # Syntax
///
/// ```rust
/// use rastair2_vcf::{format_field, FieldNumber};
/// type Type = u32; // or any other type that implements FormatFieldValue
///
/// format_field!(Name(Type), "ID", "Description", FieldNumber::OneValPerAlt);
/// ```
///
/// This will define a struct `Name` that implements the [`FormatField`] trait (as well as [`crate::VcfField`] and [`crate::HeaderField`]).
///
/// The last parameter must be a variant of [`crate::FieldNumber`].
#[macro_export]
macro_rules! format_field {
    ($name:ident($type:ty), $id:expr, $desc:expr, $number:expr) => {
        #[doc = $desc]
        #[doc = "format field for VCF output"]
        #[derive(Debug, Clone)]
        pub struct $name(pub Vec<$type>);

        impl $crate::VcfField for $name {
            const ID: &'static str = $id;
            const NUMBER: $crate::FieldNumber = $number;
        }

        impl $crate::HeaderField for $name {
            const DESCRIPTION: &'static str = $desc;
        }

        impl $crate::FormatField for $name {
            type Type = $type;

            fn write(&self, record: &mut rust_htslib::bcf::Record) -> color_eyre::Result<()> {
                use $crate::VcfField as _;

                <$type as $crate::FormatFieldValue>::write(record, Self::ID, &self.0)
            }
        }
    };
}

#[cfg(test)]
mod tests {
    use super::super::FieldNumber;
    use super::*;
    use insta::assert_snapshot;
    use rust_htslib::bcf::{Format, Header, Writer};
    use tempfile::TempDir;

    #[test]
    fn format_header() {
        format_field!(Foo(String), "GT", "Foo", FieldNumber::OnePerGenotype);

        assert_snapshot!(
            Foo::header(),
            @r###"##FORMAT=<ID=GT,Number=G,Type=String,Description="Foo">"###
        );
    }

    #[test]
    fn integers() -> Result<()> {
        format_field!(FieldU32(u32), "U32", "Test u32", FieldNumber::OneValPerAlt);
        format_field!(FieldU64(u64), "U64", "Test u64", FieldNumber::OneValPerAlt);
        format_field!(FieldI32(i32), "I32", "Test i32", FieldNumber::OneValPerAlt);
        format_field!(FieldI64(i64), "I64", "Test i64", FieldNumber::OneValPerAlt);
        format_field!(FieldUsize(usize), "Usize", "Test usize", FieldNumber::OneValPerAlt);

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
        FieldU32(vec![1]).write(&mut record)?;
        assert_eq!((*record.format(b"U32").integer()?)[0], &[1i32]);

        FieldI32(vec![42]).write(&mut record)?;
        assert_eq!((*record.format(b"I32").integer()?)[0], &[42i32]);

        FieldU64(vec![42]).write(&mut record)?;
        assert_eq!((*record.format(b"U64").integer()?)[0], &[42i32]);

        FieldI64(vec![42]).write(&mut record)?;
        assert_eq!((*record.format(b"I64").integer()?)[0], &[42i32]);

        FieldUsize(vec![42]).write(&mut record)?;
        assert_eq!((*record.format(b"Usize").integer()?)[0], &[42i32]);

        // lists
        FieldU32(vec![1, 2]).write(&mut record)?;
        assert_eq!((*record.format(b"U32").integer()?)[0], &[1i32, 2]);

        FieldI64(vec![1, 2]).write(&mut record)?;
        assert_eq!((*record.format(b"I64").integer()?)[0], &[1i32, 2]);

        vcf.write(&record).wrap_err("Failed to write record")?;
        drop(vcf);
        let result = std::fs::read_to_string(temp_file).wrap_err("read temp file")?;
        assert_snapshot!(result);

        // test that we catch overflow
        assert!(FieldU32(vec![u32::MAX]).write(&mut record).is_err());
        assert!(FieldU64(vec![u64::MAX]).write(&mut record).is_err());
        assert!(FieldI64(vec![i64::MAX]).write(&mut record).is_err());
        // but i32 is base type so it's fine
        assert!(FieldI32(vec![i32::MAX]).write(&mut record).is_ok());

        Ok(())
    }
}
