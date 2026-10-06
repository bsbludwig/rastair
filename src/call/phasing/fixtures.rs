//! Hand-built columns for the phasing tests.

use super::{ColumnPhase, FragmentId, GtAllele, PhaseCall, PhaseObservation, PhaseSet};
use crate::{
    call::variant_calling::{EstimatedGenotype, GenotypeTag},
    metrics::{Alt, AltCall, PileupMetrics, PositionMetrics},
    sequence::{ChunkRegion, Region},
    utils::default,
};
use seqair_types::Pos1;
use seqair_types::{Base, BaseQuality, Probability, Strand};
use std::{num::NonZeroU8, sync::Arc};

pub(crate) const ALT_1: NonZeroU8 = NonZeroU8::MIN;
pub(crate) const ALT_2: NonZeroU8 = NonZeroU8::new(2).expect("2 > 0");

/// A phase call in the shorthand the tests use: `set` as the VCF prints it
/// (`1001` for a block anchored on the 0-based column 1000), `0` for the
/// reference, `n` for the `n`-th alt.
pub(crate) fn call(set: u32, first: u8) -> PhaseCall {
    PhaseCall { set: PhaseSet(Pos1::new(set).expect("a 1-based phase set")), first: allele(first) }
}

pub(crate) fn allele(index: u8) -> GtAllele {
    NonZeroU8::new(index).map_or(GtAllele::Ref, GtAllele::Alt)
}

pub(crate) fn obs(fragment: u64, base: Base, strand: Strand) -> PhaseObservation {
    graded(fragment, base, strand, 30)
}

pub(crate) fn graded(fragment: u64, base: Base, strand: Strand, qual: u8) -> PhaseObservation {
    PhaseObservation {
        fragment: FragmentId::new(fragment).expect("non-zero test fragment id"),
        base,
        strand,
        qual: BaseQuality::from_byte(qual),
    }
}

/// A column with whatever genotype, alts and observations a test needs, and
/// defaults everywhere else.
pub(crate) fn column(
    pos: u32,
    reference: Base,
    alts: &[(Base, AltCall)],
    genotype: Option<GenotypeTag>,
    observations: &[PhaseObservation],
) -> PileupMetrics {
    let mut metrics = PileupMetrics {
        region: Arc::new(ChunkRegion {
            region: Region { contig: "chr_test".into(), start: 1000, end: 2000 },
            last_position: 2000,
            overlap_start: 0,
            overlap_end: 0,
        }),
        pos,
        reference_base: reference,
        context: default(),
        pos_metrics: PositionMetrics::default(),
        pos_filters: default(),
        ref_metrics: default(),
        alts: alts
            .iter()
            .map(|(base, call)| Alt {
                base: *base,
                metrics: default(),
                filters: default(),
                call: call.clone(),
            })
            .collect(),
        before_counts: default(),
        after_counts: default(),
        tags: default(),
        indel_data: None,
        phase: (!observations.is_empty())
            .then(|| Box::new(ColumnPhase::Observed(observations.into()))),
    };
    metrics.pos_metrics.extended.genotype = genotype.map(|genotype| EstimatedGenotype {
        genotype,
        likelihood: Probability::ONE,
        confidence: Probability::ONE,
    });
    metrics
}

/// The ordinary case: a `0/1` site whose single alt was called.
pub(crate) fn het(
    pos: u32,
    reference: Base,
    alt: Base,
    observations: &[PhaseObservation],
) -> PileupMetrics {
    column(
        pos,
        reference,
        &[(alt, AltCall::RealVariant)],
        Some(GenotypeTag::RefHet(ALT_1)),
        observations,
    )
}
