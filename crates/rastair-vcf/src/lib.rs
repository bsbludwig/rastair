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
//! # use rastair_types::smallvec::{SmallVec, smallvec};
//! # use std::collections::BTreeSet;
//! # use rastair_types::SmolStr;
//! # use tempfile::TempDir;
//! use rastair_vcf::*;
//!
//! info_field!(AlleleFrequency(f64), "AF", "Allele Frequency", InfoFieldNumber::OnePerAlt);
//! format_field!(ReadDepth(u32), "RD", "Read Depth", FormatFieldNumber::Num(1));
//! filter!(q10, "Quality below 10");
//! filter!(s50, "Less than 50% of samples have data");
//!
//! vcf_record!(
//!     filters: [q10, s50],
//!     info: [AlleleFrequency default],
//!     format: [ReadDepth default],
//!     min_samples: 1
//! );
//!
//! # let temp_dir = TempDir::new()?;
//! # let temp_file = temp_dir.path().join("test.vcf");
//! let writer = VcfBuilder::new(&temp_file, VcfFormat::Vcf, Compression::Off, 1)?;
//!
//! let contigs = [Contig { name: SmolStr::new("1"), length: 1000 }];
//! let samples = [SmolStr::new("sample")];
//! let mut vcf = writer.build::<Record>(&contigs, &samples)?;
//!
//! {
//!     let data = Record {
//!         main: VcfFixedFields {
//!             chrom: "1".into(),
//!             pos: 7,
//!             id: BTreeSet::from(["rs123".into()]),
//!             r#ref: "A".into(),
//!             alt: smallvec!["C".into(), "G".into()],
//!             qual: Some(50.0),
//!         },
//!         filters: { let mut f = Filters::new(); f.add(q10.filter()); f.add(s50.filter()); f },
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
//!     vcf.add(&data)?;
//! }
//! # Ok(()) }
//! ```

#![deny(missing_docs)]
#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

pub(crate) mod fields;
pub(crate) mod filters;
pub(crate) mod fixed_fields;
pub(crate) mod record;

pub mod reflect;

pub mod standard_fields;

use color_eyre::{Result, eyre::Context};
use rastair_types::SmolStr;
pub use rust_htslib::bcf::Format as VcfFormat;
use rust_htslib::bcf::{Header, Writer};
use std::{
    collections::BTreeMap,
    marker::PhantomData,
    path::{Path, PathBuf},
};

pub use crate::fields::{
    FormatField, FormatFieldNumber, FormatFieldValue, HeaderField, InfoField, InfoFieldNumber,
    InfoFieldValue, StrandSpecificInfoField, VcfField,
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
    threads: usize,
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
    pub fn new(
        target: &Path,
        format: VcfFormat,
        compressed: Compression,
        threads: usize,
    ) -> Result<Self> {
        let header = Header::new();

        Ok(Self { target: target.to_path_buf(), format, compressed, header, threads })
    }

    /// Add a line to the VCF header.
    pub fn add_header_line(&mut self, line: impl AsRef<[u8]>) {
        self.header.push_record(line.as_ref());
    }

    /// Initialize the VCF with a given record type, contigs, and samples.
    ///
    /// This will start writing to the target file specified in the builder.
    pub fn build<R: WriteToVcf>(
        mut self,
        contigs: &[Contig],
        samples: &[SmolStr],
    ) -> Result<VcfFile<R>>
    where
        R::Config: Default,
    {
        R::write_header(&mut self.header).wrap_err("Failed to write VCF header")?;

        for contig in contigs {
            let Contig { name, length } = contig;
            self.add_header_line(format!(r#"##contig=<ID={name},length={length}>"#));
        }
        for sample in samples {
            self.header.push_sample(sample.as_bytes());
        }

        let mut vcf = Writer::from_path(
            &self.target,
            &self.header,
            matches!(self.compressed, Compression::Off),
            self.format,
        )
        .wrap_err_with(|| format!("Failed to create VCF writer for `{}`", self.target.display()))?;

        let extra_threads = self.threads.saturating_sub(1); // we are one of the threads already
        if extra_threads > 0 {
            vcf.set_threads(extra_threads).wrap_err("Failed to set threads for VCF writer")?;
        }

        let mut chromosomes = BTreeMap::new();
        for contig in contigs {
            let id = vcf
                .header()
                .name2rid(contig.name.as_bytes())
                .wrap_err_with(|| format!("Failed to add contig `{}` to header", contig.name))?;
            chromosomes.insert(contig.name.clone(), id);
        }

        Ok(VcfFile {
            chromosomes,
            samples: 0,
            record_type: PhantomData,
            writer: vcf,
            field_config: Default::default(),
        })
    }
}

/// VCF writer for a specific record type
pub struct VcfFile<R: WriteToVcf> {
    /// Maps chromosome names to their IDs in the VCF header
    pub chromosomes: BTreeMap<SmolStr, u32>,
    /// Number of samples in the VCF
    pub samples: u16,
    record_type: PhantomData<R>,
    writer: Writer,
    /// Configuration for which fields to write
    field_config: R::Config,
}

impl<R: WriteToVcf> VcfFile<R> {
    /// Set a custom field configuration
    ///
    /// This allows controlling which fields are written to the VCF output.
    pub fn with_config(mut self, config: R::Config) -> Self {
        self.field_config = config;
        self
    }

    /// Get a mutable reference to the field configuration
    pub fn config_mut(&mut self) -> &mut R::Config {
        &mut self.field_config
    }

    /// Add a record to the VCF
    pub fn add(&mut self, data: &R) -> Result<()> {
        let mut record = self.writer.empty_record();
        data.write(&mut record, &self.field_config).wrap_err("Failed to write record")?;
        self.writer.write(&record).wrap_err("Failed to write record to VCF file")?;
        Ok(())
    }
}

/// Represents a contig (chromosome) in the VCF file
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, serde::Serialize, serde::Deserialize)]
pub struct Contig {
    /// Name of the contig (chromosome)
    pub name: SmolStr,
    /// Length of the contig
    pub length: u64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{fixed_fields::VcfFixedFields, vcf_record};
    use color_eyre::{Result, eyre::Context};
    use fields::InfoFieldNumber;
    use insta::assert_snapshot;
    use rastair_types::smallvec::{SmallVec, smallvec};
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
            info: [AlleleFrequency default],
            format: [Example default],
            min_samples: 1
        );

        // Use FieldConfig to suppress unused method warning
        let config = FieldConfig::default();
        let _ = config.with_field_ids(&[], &[]);

        let temp_dir = TempDir::new()?;
        let temp_file = temp_dir.path().join("test.vcf");
        let writer = VcfBuilder::new(&temp_file, VcfFormat::Vcf, Compression::Off, 1)?;

        let contigs = [Contig { name: SmolStr::new("1"), length: 1000 }];
        let samples = [SmolStr::new("sample")];
        let mut vcf = writer.build::<Record>(&contigs, &samples)?;

        {
            let data = Record {
                main: VcfFixedFields {
                    chrom: "1".into(),
                    pos: 7,
                    id: BTreeSet::from(["rs123".into()]),
                    r#ref: "A".into(),
                    alt: SmallVec::from(["C".into(), "G".into()]),
                    qual: Some(50.0),
                },
                filters: {
                    let mut f = Filters::new();
                    f.add(q10.filter());
                    f.add(s50.filter());
                    f
                },
                info: Info { allele_frequency: AlleleFrequency(smallvec![0.5, 0.75]) },
                samples: smallvec![Format {
                    example: Example(smallvec!["0/1".into(), "1/1".into()])
                }],
            };

            vcf.add(&data)?;
        }

        drop(vcf);
        let result = std::fs::read_to_string(temp_file).wrap_err("read temp file")?;
        assert_snapshot!(result);

        Ok(())
    }
}
