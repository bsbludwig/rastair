//! Methylation β per CpG a position takes part in.
//!
//! A CpG is a `C` and a `G` next to each other on one chromosome
//! (`PileupMetrics::cpg_side_of`, from the genotypes here and next door). A
//! position's CpG on a side is the reference's (`Original`) when the
//! reference base takes that side, else a called alt's (`DeNovo`); where no
//! chromosome has both bases there is none. β counts the reads on the strand
//! that converts the side whose neighbour in the read is the CpG partner
//! (`read_counts`), which leaves out a chromosome without the partner;
//! `M5mC_Strands` and the BED `unmod` and `mod` columns count the same reads.
//! Only a chromosome reading the converted base beside the partner (`TG`
//! beside a `CG`) puts its reads among them, and β then counts half as its
//! ([`confounded`], `CPGconfounded`); the CpG's other side counts none of them. The VCF writes one slot per CpG, `Original` first, and a
//! CpG no read shows has no β ([`written`]).

use crate::{
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

    let reads = read_counts(record, side);
    let Some(total_count) = NonZeroU32::new(reads.modified + reads.unmodified) else {
        return Ok(None);
    };
    let beta = beta(&reads, confounded(record, side));

    Ok(Some((
        origin,
        CpgBeta {
            beta: Probability::new(beta).this_is_a_bug()?,
            mod_count: reads.modified,
            total_count,
        },
    )))
}

/// The CpGs the VCF (`M5mC`, `DPM5mC`, `ADM5mC`) writes for this position:
/// one per CpG [`origins`], Original first, so a value's slot says which CpG
/// it measures; `None` for a CpG no read shows, which has no β. A reference
/// CpG's position keeps its Original slot, empty, where no chromosome keeps
/// the CpG.
pub fn written(record: &PileupMetrics) -> SmallVec<(CpgOrigin, Option<CpgBeta>), 2> {
    let methylated = &record.pos_metrics.methylated;
    let mut origins = origins(record);
    if record.reference_cpg().is_some() && !origins.contains(&CpgOrigin::Original) {
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
/// no chromosome has a CpG's base here on that side.
pub(crate) fn cpg_origin(record: &PileupMetrics, side: CpgSide) -> Option<CpgOrigin> {
    // The reference base takes this side, in a reference CpG or beside a
    // called alt that makes one, on a chromosome that keeps it; else a called
    // alt does.
    if record.cpg_side_of(record.ref_base()) == Some(side) {
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

/// Whether a converted read on `side` may be a chromosome's own base rather
/// than methylation: the chromosome without the CpG reads the converted base
/// (`T` for the `C`, `A` for the `G`) with the partner beside it, so β's reads
/// count it.
pub(crate) fn confounded(record: &PileupMetrics, side: CpgSide) -> bool {
    record.chromosomes().is_some() && record.reads_beside_partner(side.mod_base(), side)
}

/// β of a side's reads. At a [`confounded`] side half the reads are taken to
/// be the other chromosome's, all converted, and only the converted reads
/// beyond that half the CpG's.
fn beta(reads: &CpgReads, confounded: bool) -> f64 {
    let (modified, unmodified) = (reads.modified.f(), reads.unmodified.f());
    let modified =
        if confounded { (modified - (modified + unmodified) / 2.).max(0.0) } else { modified };
    modified / (modified + unmodified)
}
