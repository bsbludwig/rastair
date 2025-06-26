//! VCF filters
//!
//! From spec:
//!
//! > PASS if this position has passed all filters, i.e., a call is made at this position.
//! > Otherwise, if the site has not passed all filters, a semicolon-separated list of codes for filters that fail. e.g.
//! > “q10;s50” might indicate that at this site the quality is below 10 and the number of samples with data is below
//! > 50% of the total number of samples. ‘0’ is reserved and must not be used as a filter String. If filters have
//! > not been applied, then this field must be set to the MISSING value. (String, no whitespace or semicolons
//! > permitted, duplicate values not allowed.)

use smol_str::SmolStr;

/// A filter that can be applied to VCF records
pub trait VcfFilter {
    /// The name of the filter, used in the VCF header
    const NAME: &'static str;

    /// Definition of the filter
    fn header() -> String;

    /// Returns the filter name as it should appear in the VCF record
    fn filter(&self) -> SmolStr {
        SmolStr::new_static(Self::NAME)
    }
}

/// Define a VCF filter by its name and description.
///
/// This creates a struct with the filter name and implements the [`VcfFilter`] trait.
#[macro_export]
macro_rules! filter {
    ($name:ident, $description:expr) => {
        #[doc = $description]
        #[doc = "info field for VCF output"]
        #[derive(Debug, Clone)]
        #[allow(non_camel_case_types, clippy::upper_case_acronyms)]
        pub struct $name;

        impl $crate::VcfFilter for $name {
            const NAME: &'static str = stringify!($name);

            fn header() -> String {
                format!("##FILTER=<ID={},Description=\"{}\">\n", stringify!($name), $description)
            }
        }

        impl rust_htslib::bcf::record::FilterId for $name {
            fn id_from_header(
                &self,
                header: &rust_htslib::bcf::header::HeaderView,
            ) -> rust_htslib::errors::Result<rust_htslib::bcf::header::Id> {
                header.name_to_id(stringify!($name).as_bytes())
            }

            fn is_pass(&self) -> bool {
                // originally, we would check this:
                // matches!(stringify!($name).as_bytes(), b"PASS" | b".")
                // but since we are defining only custom filters, we know they never mean "pass"
                false
            }
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use color_eyre::eyre::Context;
    use insta::assert_snapshot;
    use rust_htslib::bcf::{Format, Writer, header::Header};
    use tempfile::TempDir;

    #[test]
    fn test_filter_macro() -> color_eyre::Result<()> {
        filter!(q10, "Quality below 10");
        filter!(s50, "Less than 50% of samples have data");

        let temp_dir = TempDir::new()?;
        let temp_file = temp_dir.path().join("test.vcf");

        let mut header = Header::new();
        header.push_record(b"##fileformat=VCFv4.2");
        header.push_record(br#"##contig=<ID=1,length=10>"#);

        header.push_record(q10::header().as_bytes());
        header.push_record(s50::header().as_bytes());

        let mut vcf = Writer::from_path(&temp_file, &header, true, Format::Vcf)?;
        let mut record = vcf.empty_record();

        assert!(record.has_filter("PASS".as_bytes()));

        record.push_filter(&q10)?;
        record.push_filter(&s50)?;

        vcf.write(&record)?;

        drop(vcf);
        let result = std::fs::read_to_string(temp_file).wrap_err("read temp file")?;
        assert_snapshot!(result);

        Ok(())
    }
}
