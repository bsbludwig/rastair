use color_eyre::{Result, Section as _, eyre::WrapErr as _};
use rust_htslib::bcf::Record;
use smallvec::SmallVec;
use smol_str::SmolStr;
use std::collections::BTreeSet;

/// Fixed fields in a VCF record
///
/// See [VCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf) section 1.6.1
/// from which the descriptions of the fields are taken.
pub struct VcfFixedFields {
    /// Chromosome name
    ///
    /// NOTE: Each chromosome must be registered in the VCF header using a `##contig` line.
    ///
    /// An identifier from the reference genome or an angle-bracketed ID String
    /// (`<ID>`) pointing to a contig in the assembly file (cf. the `##assembly`
    /// line in the header). All entries for a specific CHROM must form a
    /// contiguous block within the VCF file. (String, no whitespace permitted,
    /// Required)
    pub chrom: SmolStr,
    /// Position in the chromosome (1-based)
    ///
    /// Positions are sorted numerically, in increasing order, within each
    /// reference sequence CHROM. It is permitted to have multiple records with
    /// the same POS. Telomeres are indicated by using positions 0 or N+1, where
    /// N is the length of the corresponding chromosome or contig
    pub pos: u32,
    /// Identifier: Semicolon-separated list of unique identifiers where available
    ///
    /// If this is a dbSNP variant the rs number(s) should be used. No
    /// identifier should be present in more than one data record. If there is
    /// no identifier available, then the MISSING value should be used. (String,
    /// no whitespace or semicolons permitted, duplicate values not allowed.)
    pub id: BTreeSet<SmolStr>,
    /// Reference base(s)
    ///
    /// Each base must be one of A,C,G,T,N (case insensitive). Multiple bases
    /// are permitted. The value in the POS field refers to the position of the
    /// first base in the String.
    ///
    /// For simple insertions and deletions in which either the REF or one of
    /// the ALT alleles would otherwise be null/empty, the REF and ALT Strings
    /// must include the base before the variant (which must be reflected in the
    /// POS field), unless the variant occurs at position 1 on the contig in
    /// which case it must include the base after the variant; this padding base
    /// is not required (although it is permitted) e.g. for complex
    /// substitutions or other variants where all alleles have at least one base
    /// represented in their Strings. If any of the ALT alleles is a symbolic
    /// allele (an angle-bracketed ID String `<ID>`) then the padding base is
    /// required and POS denotes the coordinate of the base preceding the
    /// polymorphism. The exception to this is the <*> symbolic allele for which
    /// the reference call interval includes the POS base. Tools processing VCF
    /// files are not required to preserve case in the REF allele Strings.
    ///
    /// If the reference sequence contains IUPAC ambiguity codes not allowed by
    /// this specification (such as R = A/G), the ambiguous reference base must
    /// be reduced to a concrete base by using the one that is first
    /// alphabetically (thus R as a reference base is converted to A in VCF.)
    pub r#ref: SmolStr,
    /// Alternative base(s)
    ///
    /// Comma-separated list of alternate non-reference alleles. These alleles
    /// do not have to be called in any of the samples. Each allele in this list
    /// must be one of:
    /// - a non-empty String of bases (A,C,G,T,N; case insensitive)
    /// - the ‘*’ symbol (allele missing due to overlapping deletion)
    /// - the MISSING value ‘.’ (no variant)
    /// - an angle-bracketed ID String (`<ID>`)
    /// - the unspecified allele `<*>` as described in Section 5.5;
    /// - or a breakend replacement string as described in Section 5.4.
    ///
    /// If there are no alternative alleles, then the MISSING value must be
    /// used. Tools processing VCF files are not required to preserve case in
    /// the allele String, except for IDs, which are case sensitive.
    ///
    /// (String; no whitespace, commas, or angle-brackets are permitted in the
    /// ID String itself)
    pub alt: SmallVec<SmolStr, 2>,
    /// Quality
    ///
    /// Phred-scaled quality score for the assertion made in ALT. i.e. `−10
    /// log10` prob (call in ALT is wrong). If ALT is ‘.’ (no variant) then this
    /// is `−10 log10` prob (variant), and if ALT is not ‘.’ this is `−10 log10`
    /// prob (no variant). If unknown, the MISSING value must be specified.
    pub qual: Option<f32>,
    // Following fields are application-specific:
    // - FILTER
    // - INFO
    // - FORMAT + Samples
}

impl VcfFixedFields {
    /// Set the fixed fields in the VCF record
    pub fn write(&self, record: &mut Record) -> Result<()> {
        // PERF: Change this to keep track of the `rid` directly in the struct
        let rid = record
            .header()
            .name2rid(self.chrom.as_bytes())
            .wrap_err_with(|| format!("Failed to find chromosome {}", self.chrom))
            .note("Program error: All chromosomes need to be registered in header `##contig`")?;
        record.set_rid(Some(rid));

        record.set_pos(i64::from(self.pos));

        record.clear_id().wrap_err("Failed to reset id column")?; // Load-bearing: Push segfaults without this
        for id in &self.id {
            record.push_id(id.as_bytes()).wrap_err("Failed to push ID")?;
        }

        // Set both ref and alt alleles
        let alleles: SmallVec<_, 6> = std::iter::once(self.r#ref.as_bytes())
            .chain(self.alt.iter().map(|alt| alt.as_bytes()))
            .collect();
        record.set_alleles(alleles.as_slice()).wrap_err("Failed to set alleles")?;

        if let Some(qual) = self.qual {
            record.set_qual(qual);
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use color_eyre::Result;
    use insta::assert_snapshot;
    use rust_htslib::bcf::{Format, Writer, header::Header};
    use std::collections::BTreeSet;
    use tempfile::TempDir;

    #[test]
    fn fixed_fields() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let temp_file = temp_dir.path().join("test.vcf");

        let mut header = Header::new();
        header.push_record(br#"##contig=<ID=1,length=10>"#);

        let vcf = Writer::from_path(&temp_file, &header, true, Format::Vcf)?;
        let fields = VcfFixedFields {
            chrom: "1".into(),
            pos: 7,
            id: BTreeSet::from(["rs123".into()]),
            r#ref: "A".into(),
            alt: SmallVec::from(["C".into(), "G".into()]),
            qual: Some(50.0),
        };

        let mut record = vcf.empty_record();
        fields.write(&mut record)?;
        let vcf = record.to_vcf_string()?;
        assert_snapshot!(vcf, @"1\t8\trs123\tA\tC,G\t50\t.\t.");

        Ok(())
    }
}
