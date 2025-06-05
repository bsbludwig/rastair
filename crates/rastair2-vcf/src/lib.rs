//! VCF (Variant Call Format) writer module
//!
//! Uses `rust_htslib` for writing VCF files, so it can also handle BCF (binary
//! VCF) files as well as compressed VCF/BCF files.
//!
//! # Examples
//!
//! ```rust
//! # fn main() -> color_eyre::Result<()> {
//! # use color_eyre::eyre::WrapErr;
//! # use smallvec::{SmallVec, smallvec};
//! # use std::collections::BTreeSet;
//! # use smol_str::SmolStr;
//! # use tempfile::TempDir;
//! use rastair2_vcf::*;
//!
//! info_field!(AlleleFrequency(f64), "AF", "Allele Frequency", InfoFieldNumber::OnePerAlt);
//! format_field!(ReadDepth(u32), "RD", "Read Depth", FormatFieldNumber::Num(1));
//! filter!(q10, "Quality below 10");
//! filter!(s50, "Less than 50% of samples have data");
//!
//! vcf_record!(
//!     filters: [q10, s50],
//!     info: [AlleleFrequency],
//!     format: [ReadDepth],
//!     min_samples: 1
//! );
//!
//! # let temp_dir = TempDir::new()?;
//! # let temp_file = temp_dir.path().join("test.vcf");
//! let writer = VcfBuilder::new(&temp_file, VcfFormat::Vcf, Compression::Off)?;
//!
//! let contigs = [SmolStr::new("1")];
//! let samples = [SmolStr::new("sample")];
//! let mut vcf = writer.build::<Record>(&contigs, &samples)?;
//!
//! {
//!     let data = Record {
//!         fixed_fields: VcfFixedFields {
//!             chrom: "1".into(),
//!             pos: 7,
//!             id: BTreeSet::from(["rs123".into()]),
//!             r#ref: "A".into(),
//!             alt: smallvec!["C".into(), "G".into()],
//!             qual: 50.0,
//!         },
//!         filters: Filters::new().add(q10).add(s50),
//!         info: Info {
//!             allele_frequency: AlleleFrequency(smallvec![0.5, 0.75]),
//!         },
//!         samples: smallvec![
//!             Format {
//!                 read_depth: ReadDepth(smallvec![4])
//!             }
//!         ],
//!     };
//!
//!     vcf.add(data)?;
//! }
//! # Ok(()) }
//! ```

#![deny(missing_docs)]

pub(crate) mod fields;
pub(crate) mod filters;
pub(crate) mod fixed_fields;
pub(crate) mod record;

pub mod standard_fields;

use color_eyre::{Result, eyre::Context};
pub use rust_htslib::bcf::Format as VcfFormat;
use rust_htslib::bcf::{Header, Writer};
use smol_str::SmolStr;
use std::{
    collections::BTreeMap,
    marker::PhantomData,
    path::{Path, PathBuf},
};

pub use crate::fields::{
    FormatField, FormatFieldNumber, FormatFieldValue, HeaderField, InfoField, InfoFieldNumber,
    InfoFieldValue, VcfField,
};
pub use filters::VcfFilter;
pub use fixed_fields::VcfFixedFields;
pub use record::WriteToVcf;

/// Builder for creating a VCF writer
pub struct VcfBuilder {
    header: Header,
    target: PathBuf,
    format: VcfFormat,
    compressed: Compression,
}

/// Compression options for the VCF writer
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compression {
    /// Use compression (e.g., bgzip)
    On,
    /// Do not use compression
    Off,
}

impl VcfBuilder {
    /// Create a new VCF builder
    ///
    /// Note: This only initializes the header. It does not create the VCF file yet.
    pub fn new(target: &Path, format: VcfFormat, compressed: Compression) -> Result<Self> {
        let header = Header::new();

        Ok(Self { target: target.to_path_buf(), format, compressed, header })
    }

    /// Initialize the VCF with a given record type, contigs, and samples.
    ///
    /// This will start writing to the target file specified in the builder.
    pub fn build<R: WriteToVcf>(
        mut self,
        contigs: &[SmolStr],
        samples: &[SmolStr],
    ) -> Result<Vcf<R>> {
        R::write_header(&mut self.header).wrap_err("Failed to write VCF header")?;

        for chrom in contigs {
            self.header.push_record(format!(r#"##contig=<ID={},length=0>"#, chrom).as_bytes());
        }
        for sample in samples {
            self.header.push_sample(sample.as_bytes());
        }

        let vcf = Writer::from_path(
            &self.target,
            &self.header,
            matches!(self.compressed, Compression::Off),
            self.format,
        )
        .wrap_err_with(|| format!("Failed to create VCF writer for `{}`", self.target.display()))?;

        let mut chromosomes = BTreeMap::new();
        for chrom in contigs {
            let id = vcf
                .header()
                .name2rid(chrom.as_bytes())
                .wrap_err_with(|| format!("Failed to add contig `{}` to header", chrom))?;
            chromosomes.insert(chrom.clone(), id);
        }

        Ok(Vcf { chromosomes, samples: 0, record_type: PhantomData, writer: vcf })
    }
}

/// VCF writer for a specific record type
pub struct Vcf<R: WriteToVcf> {
    /// Maps chromosome names to their IDs in the VCF header
    pub chromosomes: BTreeMap<SmolStr, u32>,
    /// Number of samples in the VCF
    pub samples: u16,
    record_type: PhantomData<R>,
    writer: Writer,
}

impl<R: WriteToVcf> Vcf<R> {
    /// Add a record to the VCF
    pub fn add(&mut self, data: R) -> Result<()> {
        let mut record = self.writer.empty_record();
        data.write(&mut record).wrap_err("Failed to write record")?;
        self.writer.write(&record).wrap_err("Failed to write record to VCF")?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{fixed_fields::VcfFixedFields, vcf_record};
    use color_eyre::{Result, eyre::Context};
    use fields::InfoFieldNumber;
    use insta::assert_snapshot;
    use smallvec::{SmallVec, smallvec};
    use std::collections::BTreeSet;
    use tempfile::TempDir;

    #[test]
    fn all_in() -> Result<()> {
        info_field!(AlleleFrequency(f64), "AF", "Allele Frequency", InfoFieldNumber::OnePerAlt);
        format_field!(Example(String), "Ex", "Genotype", FormatFieldNumber::Num(1));
        filter!(q10, "Quality below 10");
        filter!(s50, "Less than 50% of samples have data");

        vcf_record!(
            filters: [q10, s50],
            info: [AlleleFrequency],
            format: [Example],
            min_samples: 1
        );

        let temp_dir = TempDir::new()?;
        let temp_file = temp_dir.path().join("test.vcf");
        let writer = VcfBuilder::new(&temp_file, VcfFormat::Vcf, Compression::Off)?;

        let contigs = [SmolStr::new("1")];
        let samples = [SmolStr::new("sample")];
        let mut vcf = writer.build::<Record>(&contigs, &samples)?;

        {
            let data = Record {
                fixed_fields: VcfFixedFields {
                    chrom: "1".into(),
                    pos: 7,
                    id: BTreeSet::from(["rs123".into()]),
                    r#ref: "A".into(),
                    alt: SmallVec::from(["C".into(), "G".into()]),
                    qual: 50.0,
                },
                filters: Filters::new().add(q10).add(s50),
                info: Info { allele_frequency: AlleleFrequency(smallvec![0.5, 0.75]) },
                samples: smallvec![Format {
                    example: Example(smallvec!["0/1".into(), "1/1".into()])
                }],
            };

            vcf.add(data)?;
        }

        drop(vcf);
        let result = std::fs::read_to_string(temp_file).wrap_err("read temp file")?;
        assert_snapshot!(result);

        Ok(())
    }
}
