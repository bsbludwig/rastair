//! The calls `bam` reads: a tabix-indexed BED, or an indexed VCF/BCF.
//!
//! Both say which positions are CpGs. Only a VCF carries phase (`GT` `0|1`,
//! `PS`): the BED has one row per CpG and no row for a variant elsewhere, so
//! reads get haplotype tags ([`super::haplotag`]) from a VCF alone.
//!
//! A VCF record is classified from its `CPG` and `CPGnovo` flags
//! ([`cpg_call`]), not through the BED row `convert` makes of it. A BED row
//! describes one CpG per position, by the strand its β is read on, and a
//! position can take part in two: where a homozygous alt removes a reference
//! CpG and makes a de-novo one in its place (`CGG` to `CCG`), the row may
//! describe the lost one, and reads as an SNP. The flags keep the de-novo
//! CpG, so at such positions a VCF marks reads a BED leaves unmarked; the
//! calls are the same everywhere else.

use super::haplotag::{INDEL_FLANK, PhasedSite};
use crate::{
    bed::reader::{RastairBedReader, RastairCall},
    io::{formats::FromFileExtension as _, vcf_writer::VcfFormat},
    sequence::Region,
    utils::{
        file_helpers::{FastaReader, open_fasta},
        logging::ThisIsABug as _,
    },
};
use color_eyre::{
    Section as _,
    eyre::{Context as _, OptionExt as _, Result, eyre},
};
use rustc_hash::FxHashMap;
use seqair::vcf::{
    Allele, FilterStatus, FormatFloatOpts, FormatGt, FormatInt, Info, InfoFlag, Locus, Match,
    OnlySample, Reader, Site, VcfFile, VcfReadError,
};
use seqair_types::{Base, Pos0};
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
    Bed(Box<RastairBedReader>),
    Vcf(Box<VcfCalls>),
}

impl CallsReader {
    /// Whether the calls at `path` can carry phase, and reads so get
    /// haplotype tags: a VCF or BCF, by its name.
    pub(super) fn has_phase(path: &Path) -> bool {
        Self::vcf_format(path).is_some()
    }

    fn vcf_format(path: &Path) -> Option<VcfFormat> {
        path.file_name().and_then(|name| name.to_str()).and_then(VcfFormat::from_file_extension)
    }

    pub(super) fn open(path: &Path, fasta: &Path) -> Result<Self> {
        match Self::vcf_format(path) {
            None => RastairBedReader::new(path).map(|reader| Self::Bed(Box::new(reader))),
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
            Some(VcfFormat::VcfCompressed | VcfFormat::Bcf) => {
                let file = VcfFile::open(path)
                    .wrap_err_with(|| format!("Failed to open VCF `{}`", path.display()))?;
                if !file.has_index() {
                    return Err(eyre!("The calls VCF `{}` has no index", path.display()))
                        .with_suggestion(|| {
                            format!(
                                "`call` writes a `.csi` index next to a compressed VCF or BCF; \
                                 you can create one with:\n\
                                 \n\
                                 bcftools index {}",
                                path.display()
                            )
                        });
                }
                let keys = VcfKeys::of(&file).wrap_err_with(|| {
                    format!("Failed to read the header of `{}`", path.display())
                })?;
                let fasta = open_fasta(fasta)?;
                Ok(Self::Vcf(Box::new(VcfCalls { file, keys, reader: None, fasta })))
            }
        }
    }

    /// The calls in `region`. `contig_len` bounds the reference an indel's
    /// haplotypes are spelled with.
    pub(super) fn query(&mut self, region: &Region, contig_len: u64) -> Result<Calls> {
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
                    .into_iter()
                    .map(|call| (call.pos, call.call))
                    .collect();
                Ok(Calls { by_pos, phased: Vec::new() })
            }
            Self::Vcf(calls) => calls.query(region, contig_len),
        }
    }
}

/// What `bam` reads of a VCF record: where it is, its alleles and FILTER,
/// its INFO flags, and its one sample.
type VcfRecord = (Locus, Site, Info, OnlySample);

/// The fields `bam` reads, where the header declares them.
struct VcfKeys {
    cpg: Option<InfoFlag>,
    cpg_novo: Option<InfoFlag>,
    genotype: Option<FormatGt>,
    phase_set: Option<FormatInt>,
    beta: Option<FormatFloatOpts>,
}

impl VcfKeys {
    fn of(file: &VcfFile) -> Result<Self> {
        Ok(Self {
            cpg: file.info_key_opt("CPG")?,
            cpg_novo: file.info_key_opt("CPGnovo")?,
            genotype: file.format_key_opt("GT")?,
            phase_set: file.format_key_opt("PS")?,
            beta: file.format_key_opt("M5mC")?,
        })
    }
}

/// An indexed VCF or BCF of calls.
pub(super) struct VcfCalls {
    file: VcfFile,
    keys: VcfKeys,
    /// Opened by the first query, and pointed at each later one: that keeps
    /// its file handle and the blocks it has inflated.
    reader: Option<Reader<VcfRecord>>,
    /// An indel's haplotypes are spelled into the reference after it.
    fasta: FastaReader,
}

impl VcfCalls {
    fn query(&mut self, region: &Region, contig_len: u64) -> Result<Calls> {
        let mut calls = Calls::default();
        let pos = |at: u64| u32::try_from(at).ok().and_then(Pos0::new);
        let (Some(start), Some(end)) = (pos(region.start), pos(region.end)) else {
            return Ok(calls);
        };
        if start >= end {
            return Ok(calls);
        }
        let span = seqair::vcf::Region::new(region.contig.as_str(), start..end)
            .wrap_err_with(|| format!("Failed to name `{region}` to the VCF"))?;
        let queried = match &mut self.reader {
            Some(reader) => reader.requery(span),
            // Only the records starting in the region: one that starts
            // before it belongs to the query before.
            None => self
                .file
                .reader::<VcfRecord>()
                .matching(Match::StartsIn)
                .query(span)
                .map(|reader| self.reader = Some(reader)),
        };
        match queried {
            Ok(()) => {}
            Err(VcfReadError::UnknownContig { .. }) => {
                debug!(contig = %region.contig, "Contig not in the calls VCF");
                return Ok(calls);
            }
            Err(error) => {
                return Err(error).wrap_err_with(|| format!("Failed to fetch `{region}` from VCF"));
            }
        }
        let reader = self.reader.as_mut().ok_or_eyre("no VCF reader").this_is_a_bug()?;

        while let Some((locus, site, info, sample)) =
            reader.next_record().wrap_err("Failed to read record from VCF")?
        {
            let Some(pos) = locus.start().pos0().map(Pos0::as_u32) else { continue };

            // A position's first record with a CpG flag is its call; under
            // `--all` rejected alts follow it.
            if !calls.by_pos.contains_key(&pos)
                && let Some(call) = cpg_call(&self.keys, site, info, sample)
                    .wrap_err_with(|| format!("Failed to read the CpG call at {locus}"))?
            {
                calls.by_pos.insert(pos, call);
            }

            let fasta = &mut self.fasta;
            let reference_after = |start: u64| -> Result<Vec<u8>> {
                let stop = start.saturating_add(INDEL_FLANK).min(contig_len);
                fasta.fetch_seq(&region.contig, start, stop)
            };
            if let Some(site) = phased_site(&self.keys, pos, site, sample, reference_after)
                .wrap_err_with(|| format!("Failed to read the phased genotype at {locus}"))?
            {
                calls.phased.push(site);
            }
        }
        Ok(calls)
    }
}

/// The bases of allele `index` of a record, if it is spelled in bases.
fn bases_of(site: &Site, index: u16) -> Option<&[u8]> {
    match site.alleles().get(index)? {
        Allele::Bases(bases) => Some(bases),
        _ => None,
    }
}

/// The sample's two alleles, in `GT` order, and whether they are phased;
/// `None` unless it is diploid with both alleles called.
fn diploid(keys: &VcfKeys, sample: &OnlySample) -> Result<Option<([u16; 2], bool)>> {
    let Some(key) = &keys.genotype else { return Ok(None) };
    let Some(genotype) = sample.get(key)? else { return Ok(None) };
    let mut alleles = genotype.alleles().map(|allele| allele.allele);
    let (Some(Some(first)), Some(Some(second)), None) =
        (alleles.next(), alleles.next(), alleles.next())
    else {
        return Ok(None);
    };
    Ok(Some(([first, second], genotype.is_phased())))
}

/// The phased site a VCF record is, if its sample is a phased heterozygote
/// with a `PS` and alleles spelled in bases.
fn phased_site(
    keys: &VcfKeys,
    pos: u32,
    site: &Site,
    sample: &OnlySample,
    reference_after: impl FnOnce(u64) -> Result<Vec<u8>>,
) -> Result<Option<PhasedSite>> {
    let Some(([first, second], true)) = diploid(keys, sample)? else { return Ok(None) };
    if first == second {
        return Ok(None);
    }
    let Some(set) = keys.phase_set.as_ref().map(|key| sample.get(key)).transpose()?.flatten()
    else {
        return Ok(None);
    };
    let (Some(first), Some(second)) = (bases_of(site, first), bases_of(site, second)) else {
        return Ok(None);
    };
    PhasedSite::new(pos, set, site.alleles().reference(), [first, second], reference_after)
}

/// The CpG call a VCF record makes, if it is flagged as one.
///
/// A de-novo CpG counts where its alt passed. A reference CpG neither
/// chromosome keeps (no allele of the genotype is the reference base) is an
/// SNP, unless the alt that replaces it makes a de-novo CpG there.
fn cpg_call(
    keys: &VcfKeys,
    site: &Site,
    info: &Info,
    sample: &OnlySample,
) -> Result<Option<RastairCall>> {
    let flag = |key: &Option<InfoFlag>| -> Result<bool> {
        Ok(match key {
            Some(key) => info.get(key)?,
            None => false,
        })
    };
    let in_cpg = flag(&keys.cpg)?;
    let de_novo = flag(&keys.cpg_novo)? && matches!(site.filters().status(), FilterStatus::Pass);
    if !(in_cpg || de_novo) {
        return Ok(None);
    }

    let base_of = |index: u16| match bases_of(site, index) {
        Some(&[base]) => Some(Base::from(base)),
        _ => None,
    };
    let Some(reference) = base_of(0) else { return Ok(None) };
    // Without a genotype (`./.`) nothing says the CpG is lost.
    let [from, to] = diploid(keys, sample)?
        .and_then(|([first, second], _)| Some([base_of(first)?, base_of(second)?]))
        .unwrap_or([reference; 2]);

    let methylated = match &keys.beta {
        Some(key) => {
            sample.get(key)?.and_then(|betas| betas.flatten().next()).map(|beta| beta > 0.5)
        }
        None => None,
    };

    Ok(Some(if de_novo && !in_cpg {
        RastairCall::DeNovoCpg { base: from, methylated }
    } else if !de_novo && ![from, to].contains(&reference) {
        RastairCall::Snp { from, to }
    } else {
        RastairCall::Cpg { base: from, methylated }
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::haplotag::SiteAlleles;
    use Base::{A, C, G};

    /// What `read` makes of each VCF data line (columns separated by
    /// spaces), with the record's 0-based position.
    fn read_lines<T>(
        lines: &[&str],
        mut read: impl FnMut(&VcfKeys, u32, &Site, &Info, &OnlySample) -> Result<T>,
    ) -> Result<Vec<T>> {
        let mut vcf = String::from(
            "##fileformat=VCFv4.2\n\
             ##contig=<ID=chr1,length=1000>\n\
             ##FILTER=<ID=low_ml_score,Description=\"Rejected\">\n\
             ##INFO=<ID=CPG,Number=0,Type=Flag,Description=\"Reference CpG\">\n\
             ##INFO=<ID=CPGnovo,Number=0,Type=Flag,Description=\"De-novo CpG\">\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             ##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">\n\
             ##FORMAT=<ID=M5mC,Number=.,Type=Float,Description=\"Beta\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n",
        );
        for line in lines {
            vcf.push_str(&line.replace(' ', "\t"));
            vcf.push('\n');
        }
        let file = VcfFile::from_bytes(vcf.into_bytes(), None)?;
        let keys = VcfKeys::of(&file)?;
        let mut reader = file.reader::<VcfRecord>().open()?;
        let mut read_records = Vec::new();
        while let Some((locus, site, info, sample)) = reader.next_record()? {
            let pos = locus.start().pos0().ok_or_eyre("a record before base 1")?.as_u32();
            read_records.push(read(&keys, pos, site, info, sample)?);
        }
        Ok(read_records)
    }

    #[test]
    fn a_record_is_a_cpg_by_its_flags_and_genotype() -> Result<()> {
        let calls = read_lines(
            &[
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
            ],
            |keys, _, site, info, sample| cpg_call(keys, site, info, sample),
        )?;
        // By `Debug`: the fields are there for completeness, and nothing
        // reads them.
        let calls: Vec<String> = calls.iter().map(|call| format!("{call:?}")).collect();
        assert_eq!(
            calls,
            [
                "Some(Cpg { base: C, methylated: Some(true) })",
                "Some(Cpg { base: C, methylated: Some(false) })",
                "Some(Snp { from: T, to: T })",
                "Some(Cpg { base: C, methylated: Some(true) })",
                "Some(DeNovoCpg { base: T, methylated: Some(true) })",
                "None",
                "None",
                "Some(Cpg { base: G, methylated: None })",
                "None",
            ]
        );
        Ok(())
    }

    /// The phased sites of VCF data lines, where the reference from any
    /// position on is `flank`.
    fn sites_of(lines: &[&str], flank: &[u8]) -> Result<Vec<Option<PhasedSite>>> {
        read_lines(lines, |keys, pos, site, _, sample| {
            phased_site(keys, pos, site, sample, |_| Ok(flank.to_vec()))
        })
    }

    fn snv(pos: u32, set: i32, alleles: [Base; 2]) -> Option<PhasedSite> {
        Some(PhasedSite { pos, set, alleles: SiteAlleles::Snv(alleles) })
    }

    #[test]
    fn a_phased_het_with_a_phase_set_is_a_site() -> Result<()> {
        let sites = sites_of(
            &[
                "chr1 10 . A G . PASS . GT:PS 0|1:10",
                "chr1 20 . A G . PASS . GT:PS 1|0:10",
                "chr1 30 . A C,G . PASS . GT:PS 2|1:10",
                "chr1 40 . A G . PASS . GT:PS 0/1:.",
                "chr1 50 . A G . PASS . GT:PS 1|1:10",
                "chr1 60 . A G . PASS . GT:PS 0|1:.",
                "chr1 70 . A G . PASS . GT 0|1",
                "chr1 80 . A * . PASS . GT:PS 0|1:10",
                "chr1 90 . A . . PASS . GT:PS 0/0:.",
                "chr1 95 . A G . PASS . GT:PS .|1:10",
            ],
            b"",
        )?;
        assert_eq!(
            sites,
            [
                snv(9, 10, [A, G]),
                snv(19, 10, [G, A]),
                snv(29, 10, [G, C]),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            ]
        );
        Ok(())
    }

    /// An indel site's haplotypes are spelled past the anchor: the allele's
    /// own bases, then the reference after `REF`.
    #[test]
    fn an_indel_site_spells_both_haplotypes_into_the_reference() -> Result<()> {
        let sites = sites_of(
            &[
                "chr1 10 . A ATT . PASS . GT:PS 0|1:10",
                "chr1 20 . ACC A . PASS . GT:PS 1|0:10",
                "chr1 30 . ACC A,AGGCC . PASS . GT:PS 1|2:10",
                "chr1 40 . AC GT . PASS . GT:PS 0|1:10",
            ],
            b"gact",
        )?;
        let indel = |pos, haplotypes: [&[u8]; 2]| {
            Some(PhasedSite {
                pos,
                set: 10,
                alleles: SiteAlleles::Indel(haplotypes.map(<[u8]>::to_vec)),
            })
        };
        assert_eq!(
            sites,
            [
                indel(9, [b"GACT", b"TTGACT"]),
                indel(19, [b"GACT", b"CCGACT"]),
                indel(29, [b"GACT", b"GGCCGACT"]),
                // No shared anchor base to compare from.
                None,
            ]
        );
        Ok(())
    }
}
