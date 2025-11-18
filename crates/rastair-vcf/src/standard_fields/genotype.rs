//! This module provides the implementation of the `GenotypeAllele` format field
//! for VCF files. It uses the data from the `rust_htslib` crate to handle
//! genotype information in VCF records.
//!
//! There are multiple ways of writing genotype information. Rastair1 uses `CC`,
//! `CT`, and `TT` to show which alleles are present. In VCF, this would be
//! written `0/0`, `0/1`, and `1/1` respectively. Htslib uses
//! [`GenotypeAllele`], which also differentiates between phased and unphased
//! genotypes (i.e. knowing which parent the allele came from). Here's a quick
//! table of the mapping, assuming unphased:
//!
//! | Rastair1 | VCF | Htslib GenotypeAllele        |
//! | -------- | --- | ---------------------------- |
//! | CC       | 0/0 | `[Unphased(0)]               |
//! | CT       | 0/1 | `[Unphased(0), Unphased(1)]` |
//! | TT       | 1/1 | `[Unphased(1)]`              |
//! | TC       | 1/0 | `[Unphased(1), Unphased(0)]` |

use crate::{FormatField, FormatFieldNumber, FormatFieldValue, HeaderField, VcfField};
use color_eyre::eyre::{Context, Result};
use cstr8::CStr8;
use rastair_types::SmallVec;
pub use rust_htslib::bcf::record::GenotypeAllele as HtslibGenotypeAllele;

/// Represents a genotype in a VCF record, which can contain multiple alleles.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct Genotype(pub SmallVec<GenotypeAllele, 4>);

/// Represents a single allele in a genotype, which can be phased or unphased.
// Copy from `rust_htslib::bcf::record::GenotypeAllele` for adding serde derives
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub enum GenotypeAllele {
    /// Unphased allele with index `i`
    Unphased(i32),
    /// Phased allele with index `i`
    Phased(i32),
    /// Unphased missing allele
    UnphasedMissing,
    /// Phased missing allele
    PhasedMissing,
}

impl Genotype {
    /// Checks if the genotype is heterozygous (i.e., contains different alleles).
    pub fn heterozygous(&self) -> bool {
        match self.0.as_slice() {
            [GenotypeAllele::Unphased(i), GenotypeAllele::Unphased(j)]
            | [GenotypeAllele::Phased(i), GenotypeAllele::Phased(j)] => i != j,
            _ => false,
        }
    }

    /// Checks if the genotype is homozygous (i.e., contains the same alleles).
    pub fn homozygous(&self) -> bool {
        !self.heterozygous()
    }

    /// Checks if the genotype is homozygous and not the reference, i.e. this is SNP
    pub fn homozygous_not_ref(&self) -> bool {
        match self.0.as_slice() {
            [GenotypeAllele::Unphased(i), GenotypeAllele::Unphased(j)]
            | [GenotypeAllele::Phased(i), GenotypeAllele::Phased(j)] => i == j && *i != 0,
            _ => false,
        }
    }
}

impl From<GenotypeAllele> for HtslibGenotypeAllele {
    fn from(allele: GenotypeAllele) -> Self {
        match allele {
            GenotypeAllele::Unphased(i) => HtslibGenotypeAllele::Unphased(i),
            GenotypeAllele::Phased(i) => HtslibGenotypeAllele::Phased(i),
            GenotypeAllele::UnphasedMissing => HtslibGenotypeAllele::UnphasedMissing,
            GenotypeAllele::PhasedMissing => HtslibGenotypeAllele::PhasedMissing,
        }
    }
}

impl From<HtslibGenotypeAllele> for GenotypeAllele {
    fn from(allele: HtslibGenotypeAllele) -> Self {
        match allele {
            HtslibGenotypeAllele::Unphased(i) => GenotypeAllele::Unphased(i),
            HtslibGenotypeAllele::Phased(i) => GenotypeAllele::Phased(i),
            HtslibGenotypeAllele::UnphasedMissing => GenotypeAllele::UnphasedMissing,
            HtslibGenotypeAllele::PhasedMissing => GenotypeAllele::PhasedMissing,
        }
    }
}

impl VcfField for Genotype {
    const ID: &'static cstr8::CStr8 = cstr8::cstr8!("GT");
}

impl HeaderField for Genotype {
    const DESCRIPTION: &'static str = "Genotype";
}

impl FormatField for Genotype {
    type Type = GenotypeAllele;
    /// Number is one because the format uses specical syntax in VCF
    const NUMBER: FormatFieldNumber = FormatFieldNumber::Num(1);

    fn write(&self, record: &mut rust_htslib::bcf::Record) -> Result<()> {
        <Self::Type as FormatFieldValue>::write(record, Self::ID, self.0.as_slice())
    }
}

impl FormatFieldValue for GenotypeAllele {
    const TYPE_NAME: &'static str = "String";

    fn write(record: &mut rust_htslib::bcf::Record, _tag: &CStr8, values: &[Self]) -> Result<()> {
        record
            .push_genotypes(
                &values
                    .iter()
                    .map(|x| HtslibGenotypeAllele::from(x.clone()))
                    .collect::<SmallVec<_, 6>>(),
            )
            .wrap_err("Failed to set genotype")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use insta::assert_snapshot;
    use rastair_types::smallvec::smallvec;
    use rust_htslib::bcf::{Format, Header, Writer};
    use tempfile::TempDir;

    #[test]
    fn basic_genotype_info() -> Result<()> {
        let mut header = Header::new();
        header.push_record(b"##fileformat=VCFv4.2");
        header.push_record(br#"##contig=<ID=1,length=10>"#);
        header.push_record(Genotype::header().as_bytes());
        header.push_sample(b"one");

        let temp_dir = TempDir::new()?;
        let temp_file = temp_dir.path().join("test.vcf");
        let mut vcf = Writer::from_path(&temp_file, &header, true, Format::Vcf)?;

        {
            let mut record = vcf.empty_record();
            let genotype =
                Genotype(smallvec![GenotypeAllele::Phased(0), GenotypeAllele::Phased(1)]);
            genotype.write(&mut record)?;
            vcf.write(&record)?;
        }
        {
            let mut record = vcf.empty_record();
            let genotype =
                Genotype(smallvec![GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(0)]);
            genotype.write(&mut record)?;
            vcf.write(&record)?;
        }

        drop(vcf);
        let result = std::fs::read_to_string(temp_file).wrap_err("read temp file")?;
        assert_snapshot!(result);
        Ok(())
    }
}
