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

use color_eyre::{Result, eyre::Context as _};
use rastair_types::Base;
use smallvec::SmallVec;
use smol_str::SmolStr;

/// A list of filters
pub type FilterList = SmallVec<SmolStr, 2>;

/// Container for VCF filters holding both filters for entire records and
/// per-allele filters.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct FilterContainer {
    /// Filters that apply to the whole record
    all: FilterList,
    /// Filters that apply per allele
    per_allele: SmallVec<(Base, FilterList), 2>,
}

impl FilterContainer {
    /// Creates a new, empty filter container
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a filter that applies to the whole record
    pub fn add(&mut self, filter: SmolStr) {
        self.all.push(filter);
    }

    /// Adds a filter that applies to the whole record
    pub fn add_all(&mut self, filter: impl VcfFilter) {
        self.all.push(filter.filter());
    }

    /// Adds a filter for a specific allele
    pub fn add_per_allele(&mut self, allele: Base, filter: impl VcfFilter) {
        if let Some((_, filters)) = self.per_allele.iter_mut().find(|(a, _)| *a == allele) {
            filters.push(filter.filter());
        } else {
            self.per_allele.push((allele, {
                let mut filters = FilterList::new();
                filters.push(filter.filter());
                filters
            }));
        }
    }

    /// Checks if the given allele passes all filters
    pub fn pass_alt(&self, allele: Base) -> bool {
        let all_pass = self.all.is_empty()
            || (self.all.len() == 1 && self.all.first().expect("1 filter exists") == "PASS");
        let allele_pass =
            match self.per_allele.iter().find(|(a, _)| *a == allele).map(|(_, f)| f.as_slice()) {
                Some([]) => true,
                Some([filter]) => filter == "PASS",
                None => true,
                _ => false,
            };
        all_pass && allele_pass
    }

    /// Checks if all record-level filters pass
    fn pass_all(&self) -> bool {
        self.all.is_empty()
            || (self.all.len() == 1 && self.all.first().expect("1 filter exists") == "PASS")
    }

    /// Checks if all alleles pass all filters
    pub fn pass(&self) -> bool {
        self.pass_all() && self.per_allele.iter().all(|(allele, _)| self.pass_alt(*allele))
    }

    /// Clears all filters
    pub fn clear(&mut self) {
        self.all.clear();
        self.per_allele.clear();
    }

    /// Write to BCF record
    pub fn write_to_record(&self, record: &mut rust_htslib::bcf::Record) -> Result<()> {
        if self.all.is_empty() && self.per_allele.is_empty() {
            record.set_filters::<[u8]>(&[]).wrap_err("Failed to clear filters")?;
        } else {
            self.all
                .iter()
                .chain(self.per_allele.iter().flat_map(|(_allele, filters)| filters))
                .try_for_each(|filter| {
                    record
                        .push_filter(filter.as_bytes())
                        .wrap_err_with(|| format!("Failed to push filter {filter}"))
                })?;
        }

        Ok(())
    }
}

/// A filter that can be applied to VCF records
pub trait VcfFilter: Default {
    /// The name of the filter, used in the VCF header
    const NAME: &'static str;
    /// Description of the filter, used in the VCF header
    const DESCRIPTION: &'static str;

    /// Definition of the filter
    fn header() -> String;

    /// Returns the filter name as it should appear in the VCF record
    fn filter(&self) -> SmolStr {
        SmolStr::new_static(Self::NAME)
    }

    /// Returns the filter description
    fn description() -> crate::reflect::Filter {
        crate::reflect::Filter {
            name: SmolStr::new_static(Self::NAME),
            description: SmolStr::new_static(Self::DESCRIPTION),
        }
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
        #[derive(Debug, Clone, Default)]
        #[allow(non_camel_case_types, clippy::upper_case_acronyms)]
        pub struct $name;

        impl $name {
            const IS_PASS: bool = matches!(stringify!($name).as_bytes(), b"PASS");
        }

        impl $crate::VcfFilter for $name {
            const NAME: &'static str = stringify!($name);
            const DESCRIPTION: &'static str = $description;

            fn header() -> String {
                format!("##FILTER=<ID={},Description=\"{}\">\n", Self::NAME, Self::DESCRIPTION)
            }
        }

        impl rust_htslib::bcf::record::FilterId for $name {
            fn id_from_header(
                &self,
                header: &rust_htslib::bcf::header::HeaderView,
            ) -> rust_htslib::errors::Result<rust_htslib::bcf::header::Id> {
                header.name_to_id(cstr8::cstr8!(stringify!($name)))
            }

            fn is_pass(&self) -> bool {
                Self::IS_PASS
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
