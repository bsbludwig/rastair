//! Methylation β per CpG a position takes part in.
//!
//! A position's CpG on a side is the reference's (`Original`) when the
//! reference base takes that side and a chromosome keeps it, else a called
//! alt's (`DeNovo`); a CpG whose C and G no chromosome has together is none
//! (`process::mark_called_cpgs`). β counts the reads on the strand that
//! converts the side whose neighbour in the read is the CpG partner
//! (`read_counts`); `M5mC_Strands` and the BED `unmod` and `mod` columns count
//! the same reads. A het genotype one of whose alleles is the methylated
//! read's base, on a chromosome with the partner beside it
//! (`PileupMetrics::unconfounded`), adjusts β. The VCF writes one slot per
//! CpG, `Original` first, and a CpG no read shows has no β ([`written`]); a
//! reference CpG's position keeps its slot, empty, when it is none.

use crate::{
    call::variant_calling::GenotypeTag,
    metrics::{PileupMetrics, ReadKey},
    utils::{IntoF64, logging::ThisIsABug},
    vcf::{CpgBeta, CpgOrigin, CpgSide, Methylated},
};
use color_eyre::{Result, eyre::Context};
use seqair_types::{Probability, SmallVec};
use std::num::NonZeroU32;
use tracing::instrument;

#[cfg(test)]
mod tests;

#[instrument(
    level="debug",
    skip_all,
    fields(contig = %current.contig_name(), pos = %current.pos()),
    name = "methylation_call"
)]
/// The β of each CpG side with reads; of two de-novo sides (a called `C` and a
/// called `G` alt), the C side's.
pub fn call(current: &PileupMetrics) -> Result<Methylated> {
    let mut methylated = Methylated::default();
    for side in [CpgSide::C, CpgSide::G] {
        if let Some((origin, beta)) =
            compute_beta(current, side).wrap_err_with(|| format!("{side:?}-side beta"))?
        {
            methylated.get_mut(origin).get_or_insert(beta);
        }
    }
    Ok(methylated)
}

/// The CpG this position takes on `side`, and its methylation beta if a read
/// shows it.
fn compute_beta(record: &PileupMetrics, side: CpgSide) -> Result<Option<(CpgOrigin, CpgBeta)>> {
    let Some(origin) = cpg_origin(record, side) else { return Ok(None) };

    let CpgReads { modified: raw_mod, unmodified: raw_unmod } = read_counts(record, side);
    let Some(total_count) = NonZeroU32::new(raw_mod + raw_unmod) else { return Ok(None) };

    let adjustment = genotype_adjustment(record, side);
    let beta = adjusted_beta(raw_mod.f(), raw_unmod.f(), adjustment);

    Ok(Some((
        origin,
        CpgBeta { beta: Probability::new(beta).this_is_a_bug()?, mod_count: raw_mod, total_count },
    )))
}

/// The CpGs the VCF (`M5mC`, `DPM5mC`, `ADM5mC`) writes for this position:
/// one per CpG [`origins`], Original first, so a value's slot says which CpG
/// it measures; `None` for a CpG no read shows, which has no β. A reference
/// CpG's position keeps its Original slot, empty, where called variants leave
/// no chromosome with the CpG.
pub fn written(record: &PileupMetrics) -> SmallVec<(CpgOrigin, Option<CpgBeta>), 2> {
    let methylated = &record.pos_metrics.methylated;
    let mut origins = origins(record);
    if record.lost_reference_cpg() {
        origins.insert(0, CpgOrigin::Original);
    }
    origins.into_iter().map(|origin| (origin, methylated.get(origin).copied())).collect()
}

/// The CpGs this position takes part in, Original first.
pub fn origins(record: &PileupMetrics) -> SmallVec<CpgOrigin, 2> {
    let present: SmallVec<CpgOrigin, 2> =
        [CpgSide::C, CpgSide::G].into_iter().filter_map(|side| cpg_origin(record, side)).collect();
    [CpgOrigin::Original, CpgOrigin::DeNovo]
        .into_iter()
        .filter(|origin| present.contains(origin))
        .collect()
}

/// The CpG whose reads are this column's methylation evidence, and the side
/// it takes: the reference base's (`Original`) over a de-novo one, unless
/// only the de-novo side has reads.
pub fn evidence(record: &PileupMetrics) -> Option<(CpgSide, CpgOrigin)> {
    let has_reads = |side| {
        let CpgReads { modified, unmodified } = read_counts(record, side);
        modified + unmodified > 0
    };
    [CpgSide::C, CpgSide::G]
        .into_iter()
        .filter_map(|side| Some((side, cpg_origin(record, side)?)))
        .min_by_key(|&(side, origin)| (!has_reads(side), origin != CpgOrigin::Original))
}

/// Determine the CpG origin for this position on the given side, or None if
/// this position doesn't have a CpG allele on that side.
fn cpg_origin(record: &PileupMetrics, side: CpgSide) -> Option<CpgOrigin> {
    // The reference base takes this side, in a reference CpG or beside a
    // called alt that makes one, where a chromosome keeps it (not under a
    // `1/1`, or a `1/2` of two alts); else a called alt does.
    let keeps_reference = record
        .pos_metrics
        .genotype
        .as_ref()
        .is_none_or(|gt| matches!(gt.genotype, GenotypeTag::HomRef | GenotypeTag::RefHet(_)));
    if keeps_reference && record.cpg_side_of(record.ref_base()) == Some(side) {
        Some(CpgOrigin::Original)
    } else {
        let mut called = record.alts.iter().filter(|alt| alt.called());
        called.any(|alt| record.cpg_side_of(alt.base) == Some(side)).then_some(CpgOrigin::DeNovo)
    }
}

/// A CpG side's reads that show it methylated and unmethylated.
pub(crate) struct CpgReads {
    pub(crate) modified: u32,
    pub(crate) unmodified: u32,
}

/// The reads β counts on `side`: on the converting strand, those whose
/// neighbour in the read is the CpG partner.
pub(crate) fn read_counts(record: &PileupMetrics, side: CpgSide) -> CpgReads {
    let counts = match side {
        CpgSide::C => &record.after_counts,
        CpgSide::G => &record.before_counts,
    };
    let (strand, adj) = (side.strand(), side.adjacent_base());

    CpgReads {
        modified: counts.get(ReadKey { strand, current: side.mod_base(), adj }),
        unmodified: counts.get(ReadKey { strand, current: side.unmod_base(), adj }),
    }
}

/// How the genotype at this position affects the beta calculation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GenotypeAdjustment {
    /// Normal beta: mod / (mod + unmod).
    None,
    /// The confounding base (T for C-side, A for G-side) is present as a het
    /// allele, so some mod-base reads come from the SNP, not methylation.
    HetConfounded,
}

fn genotype_adjustment(record: &PileupMetrics, side: CpgSide) -> GenotypeAdjustment {
    let Some(gt) = &record.pos_metrics.genotype else {
        return GenotypeAdjustment::None;
    };

    // A het one of whose alleles is the base a methylated read shows: the
    // reference (`T` under a de-novo `C`) only where a chromosome carries it,
    // and only where that chromosome has the partner beside it, as β counts
    // only the reads that do.
    if gt.genotype.is_heterozygous()
        && !record.unconfounded.contains(side)
        && gt
            .genotype
            .alleles()
            .into_iter()
            .any(|allele| record.allele_base(allele) == Some(side.mod_base()))
    {
        return GenotypeAdjustment::HetConfounded;
    }

    GenotypeAdjustment::None
}

fn adjusted_beta(mod_count: f64, unmod_count: f64, adjustment: GenotypeAdjustment) -> f64 {
    match adjustment {
        GenotypeAdjustment::HetConfounded => {
            let total = mod_count + unmod_count;
            let excess_mod = (mod_count - total / 2.).max(0.0);
            excess_mod / (unmod_count + excess_mod)
        }
        GenotypeAdjustment::None => mod_count / (mod_count + unmod_count),
    }
}
