//! The calls `bam` reads: a tabix-indexed BED, or an indexed VCF/BCF.
//!
//! Both say which positions are CpGs. Only a VCF carries phase (`GT` `0|1`,
//! `PS`): the BED has one row per CpG and no row for a variant elsewhere, so
//! reads get haplotype tags ([`super::haplotag`]) from a VCF alone.
//!
//! A VCF record is classified from its `CPG` and `CPGnovo` flags
//! ([`cpg_call`]), not through the BED row `convert` makes of it: that row's
//! strand column is inferred from the record alone and differs from
//! `call --bed`'s where a position is both a lost reference CpG and a de-novo
//! one, which turns the de-novo CpG into an SNP. The flags give the calls of
//! the `call --bed` file.

use super::haplotag::{INDEL_FLANK, PhasedSite};
use crate::{
    bed::reader::{RastairBedReader, RastairCall},
    io::{formats::FromFileExtension as _, vcf_writer::VcfFormat},
    sequence::Region,
    utils::{file_helpers::FastaReader, logging::ThisIsABug as _},
};
use color_eyre::{
    Section as _,
    eyre::{Context as _, Result, eyre},
};
use rust_htslib::bcf::{
    self, Read as _,
    record::{GenotypeAllele, Numeric as _},
};
use rustc_hash::FxHashMap;
use seqair_types::Base;
use std::path::Path;
use tracing::debug;

/// What the calls file says about one stretch of a contig.
#[derive(Debug, Default)]
pub(super) struct Calls {
    /// By 0-based position.
    pub by_pos: FxHashMap<u32, RastairCall>,
    /// Sorted by position; empty for a BED.
    pub phased: Vec<PhasedSite>,
}

pub(super) enum CallsReader {
    Bed(RastairBedReader),
    Vcf(bcf::IndexedReader),
}

impl CallsReader {
    pub(super) fn open(path: &Path) -> Result<Self> {
        let format = path
            .file_name()
            .and_then(|name| name.to_str())
            .and_then(VcfFormat::from_file_extension);
        match format {
            None => RastairBedReader::new(path).map(Self::Bed),
            Some(VcfFormat::Vcf) => Err(eyre!(
                "Calls in an uncompressed VCF cannot be read by region: `{}`",
                path.display()
            ))
            .with_suggestion(|| {
                format!(
                    "Compress and index it with:\n\
                     \n\
                     bgzip {path}\n\
                     bcftools index {path}.gz",
                    path = path.display()
                )
            }),
            Some(VcfFormat::VcfCompressed | VcfFormat::Bcf) => bcf::IndexedReader::from_path(path)
                .wrap_err_with(|| format!("Failed to open indexed VCF `{}`", path.display()))
                .with_suggestion(|| {
                    format!(
                        "`call` writes a `.csi` index next to a compressed VCF or BCF; \
                         you can create one with:\n\
                         \n\
                         bcftools index {}",
                        path.display()
                    )
                })
                .map(Self::Vcf),
        }
    }

    /// Whether the calls can carry phase, and reads so get haplotype tags.
    pub(super) const fn has_phase(&self) -> bool {
        matches!(self, Self::Vcf(_))
    }

    /// The calls in `region`. `contig_len` bounds the reference an indel's
    /// haplotypes are spelled with.
    pub(super) fn query(
        &mut self,
        region: &Region,
        fasta: &FastaReader,
        contig_len: u64,
    ) -> Result<Calls> {
        match self {
            Self::Bed(reader) => {
                let region = region
                    .clone()
                    .try_into()
                    .wrap_err("Failed to convert region representation")
                    .this_is_a_bug()?;
                let by_pos = reader
                    .query(&region)
                    .wrap_err("failed to query calls file")?
                    .iter()
                    .map(|call| (call.pos, call.call.clone()))
                    .collect();
                Ok(Calls { by_pos, phased: Vec::new() })
            }
            Self::Vcf(reader) => query_vcf(reader, region, fasta, contig_len),
        }
    }
}

fn query_vcf(
    reader: &mut bcf::IndexedReader,
    region: &Region,
    fasta: &FastaReader,
    contig_len: u64,
) -> Result<Calls> {
    let mut calls = Calls::default();
    let Ok(rid) = reader.header().name2rid(region.contig.as_bytes()) else {
        debug!(contig = %region.contig, "Contig not in the calls VCF");
        return Ok(calls);
    };
    let Some(last) = region.end.checked_sub(1).filter(|last| *last >= region.start) else {
        return Ok(calls);
    };
    match reader.fetch(rid, region.start, Some(last)) {
        Ok(()) => {}
        // The index lists only contigs that hold records, so a declared
        // contig without any cannot be sought.
        Err(bcf::BcfError::GenomicSeek { .. }) => return Ok(calls),
        Err(error) => {
            return Err(error).wrap_err_with(|| format!("Failed to fetch `{region}` from VCF"));
        }
    }

    let mut record = reader.empty_record();
    while let Some(result) = reader.read(&mut record) {
        result.wrap_err("Failed to read record from VCF")?;
        // htslib also returns a record that starts before the region and
        // reaches into it.
        let Some(pos) =
            u32::try_from(record.pos()).ok().filter(|&pos| u64::from(pos) >= region.start)
        else {
            continue;
        };

        // A position's first record with a CpG flag is its call; under
        // `--all` rejected alts follow it.
        if !calls.by_pos.contains_key(&pos)
            && let Some(call) = cpg_call(&record)
        {
            calls.by_pos.insert(pos, call);
        }

        let reference_after = |start: u64| -> Result<Vec<u8>> {
            let stop = start.saturating_add(INDEL_FLANK).min(contig_len);
            if start >= stop {
                return Ok(Vec::new());
            }
            fasta.fetch_seq(&region.contig, start, stop)
        };
        if let Some(site) = PhasedSite::from_record(&record, reference_after)
            .wrap_err("Failed to read a phased genotype from a VCF record")?
        {
            calls.phased.push(site);
        }
    }
    Ok(calls)
}

/// The CpG call a VCF record makes, if it is flagged as one.
///
/// A de-novo CpG counts where its alt passed. A reference CpG neither
/// chromosome keeps (no allele of the genotype is the reference base) is an
/// SNP, unless the alt that replaces it makes a de-novo CpG there.
fn cpg_call(record: &bcf::Record) -> Option<RastairCall> {
    let in_cpg = record.info(b"CPG").flag().unwrap_or(false);
    let de_novo =
        record.info(b"CPGnovo").flag().unwrap_or(false) && record.has_filter(b"PASS".as_slice());
    if !(in_cpg || de_novo) {
        return None;
    }

    let alleles = record.alleles();
    let base_of = |index: usize| match alleles.get(index) {
        Some(&[base]) => Some(Base::from(base)),
        _ => None,
    };
    let reference = base_of(0)?;
    // Without a genotype (`./.`) nothing says the CpG is lost.
    let [from, to] = record
        .genotypes()
        .ok()
        .and_then(|genotypes| {
            let genotype = genotypes.get(0);
            let allele = |allele: &GenotypeAllele| base_of(usize::try_from(allele.index()?).ok()?);
            match &genotype[..] {
                [first, second] => Some([allele(first)?, allele(second)?]),
                _ => None,
            }
        })
        .unwrap_or([reference; 2]);

    let methylated = record.format(b"M5mC").float().ok().and_then(|betas| {
        betas.first()?.iter().find(|beta| !beta.is_missing()).map(|beta| *beta > 0.5)
    });

    Some(if de_novo && !in_cpg {
        RastairCall::DeNovoCpg { base: from, methylated }
    } else if !de_novo && ![from, to].contains(&reference) {
        RastairCall::Snp { from, to }
    } else {
        RastairCall::Cpg { base: from, methylated }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use Base::{C, G, T};

    /// The CpG calls of VCF data lines, columns separated by spaces.
    fn calls_of(lines: &[&str]) -> Result<Vec<Option<RastairCall>>> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("calls.vcf");
        let mut vcf = String::from(
            "##fileformat=VCFv4.2\n\
             ##contig=<ID=chr1,length=1000>\n\
             ##FILTER=<ID=low_ml_score,Description=\"Rejected\">\n\
             ##INFO=<ID=CPG,Number=0,Type=Flag,Description=\"Reference CpG\">\n\
             ##INFO=<ID=CPGnovo,Number=0,Type=Flag,Description=\"De-novo CpG\">\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             ##FORMAT=<ID=M5mC,Number=.,Type=Float,Description=\"Beta\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n",
        );
        for line in lines {
            vcf.push_str(&line.replace(' ', "\t"));
            vcf.push('\n');
        }
        std::fs::write(&path, vcf)?;
        let mut reader = bcf::Reader::from_path(&path)?;
        reader.records().map(|record| Ok(cpg_call(&record?))).collect()
    }

    #[test]
    fn a_record_is_a_cpg_by_its_flags_and_genotype() -> Result<()> {
        let calls = calls_of(&[
            "chr1 10 . C . . PASS CPG GT:M5mC 0/0:0.9",
            "chr1 20 . C T . PASS CPG GT:M5mC 0/1:0.2",
            // Neither chromosome keeps the reference CpG.
            "chr1 30 . C T . PASS CPG GT:M5mC 1/1:0",
            // The alt that removes the reference CpG makes a de-novo one
            // (`CGG` to `CCG`): still a CpG position.
            "chr1 40 . G C . PASS CPG;CPGnovo GT:M5mC 1/1:.,0.9",
            "chr1 50 . T C . PASS CPGnovo GT:M5mC 0/1:0.8",
            // A rejected alt makes no CpG.
            "chr1 60 . T C . low_ml_score CPGnovo GT:M5mC 0/0:.",
            "chr1 70 . A G . PASS . GT 0/1",
            // No read, no genotype: nothing says the CpG is lost.
            "chr1 80 . G . . PASS CPG GT:M5mC ./.:.",
            // An indel record at a CpG's position is not the CpG's record.
            "chr1 90 . CG C . PASS . GT 0/1",
        ])?;
        let describe = |call: &Option<RastairCall>| match call {
            Some(RastairCall::Cpg { base, methylated }) => format!("cpg {base} {methylated:?}"),
            Some(RastairCall::DeNovoCpg { base, methylated }) => {
                format!("de-novo {base} {methylated:?}")
            }
            Some(RastairCall::Snp { from, to }) => format!("snp {from}/{to}"),
            None => "none".to_owned(),
        };
        assert_eq!(
            calls.iter().map(describe).collect::<Vec<_>>(),
            [
                format!("cpg {C} Some(true)"),
                format!("cpg {C} Some(false)"),
                format!("snp {T}/{T}"),
                format!("cpg {C} Some(true)"),
                format!("de-novo {T} Some(true)"),
                "none".to_owned(),
                "none".to_owned(),
                format!("cpg {G} None"),
                "none".to_owned(),
            ]
        );
        Ok(())
    }
}
