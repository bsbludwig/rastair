//! Hand-built columns for unit tests: a column with the alts, genotype and
//! pair-HMM scores a test needs, and defaults everywhere else.

use crate::{
    call::test_helpers::pos0,
    call::{
        pair_hmm::{Allele, Cluster, Event, HaplotypeScores, ReadScores},
        pileup::indels::{IndelAllele, IndelAlleleCounts, IndelCounts, IndelData, Scored},
        variant_calling::{EstimatedGenotype, GenotypeSource, GenotypeTag},
    },
    metrics::{Alt, AltVerdict, PileupMetrics, PositionMetrics},
    sequence::{ChunkRegion, Region},
    utils::default,
};
use seqair_types::{Base, Phred, Probability, RootMeanSquare};
use std::sync::Arc;

pub(crate) fn mapq(mapq: u8) -> RootMeanSquare {
    [f64::from(mapq)].into_iter().collect()
}

/// A column with whatever genotype and alts a test needs, scored by the
/// pair-HMM over `reads` against the reference and each alt's haplotype in
/// that record order (so `scored(f, [ref, alt1, alt2])` rows), and defaults
/// everywhere else. Without reads nothing was scored.
pub(crate) fn column(
    pos: u32,
    reference: Base,
    alts: &[(Base, bool)],
    genotype: Option<GenotypeTag>,
    reads: &[ReadScores],
) -> PileupMetrics {
    let snv_scores = (!reads.is_empty() && !alts.is_empty()).then(|| {
        let events =
            alts.iter().map(|(base, _)| Event { anchor: pos0(pos), allele: Allele::Base(*base) });
        let cluster = Cluster::new(events.collect()).expect("a cluster");
        Arc::new(HaplotypeScores::new(cluster, reads.into()).expect("reads the haplotypes explain"))
    });
    let mut metrics = PileupMetrics {
        region: Arc::new(ChunkRegion {
            region: Region { contig: "chr_test".into(), start: 1000, end: 2000 },
            last_position: 2000,
            overlap_start: 0,
            overlap_end: 0,
        }),
        pos: pos0(pos),
        reference_base: reference,
        context: default(),
        pos_metrics: PositionMetrics::default(),
        beside_called_alts: default(),
        lost_to_called_alts: default(),
        ref_metrics: default(),
        alts: alts
            .iter()
            .map(|(base, call)| Alt {
                base: *base,
                metrics: default(),
                verdict: if *call {
                    AltVerdict::Called(Probability::ONE)
                } else {
                    AltVerdict::Unscored
                },
            })
            .collect(),
        before_counts: default(),
        after_counts: default(),
        indel_data: None,
        snv_scores,
        phase: None,
    };
    // Mapped well enough that no site is guarded unless a test says so.
    metrics.ref_metrics.mapq = mapq(60);
    for alt in &mut metrics.alts {
        alt.metrics.mapq = mapq(60);
    }
    metrics.pos_metrics.genotype = genotype.map(|genotype| EstimatedGenotype {
        genotype,
        quality: Phred::from_phred(0),
        source: GenotypeSource::Binomial,
    });
    metrics
}

/// Give `column` the pair-HMM's scores of `reads` against the reference and
/// the haplotypes carrying `scored`, in that record order, and the indel calls
/// it makes of them.
pub(crate) fn with_indels(
    mut column: PileupMetrics,
    scored: &[IndelAllele],
    reads: &[ReadScores],
) -> PileupMetrics {
    let events =
        scored.iter().map(|a| Event { anchor: column.pos, allele: a.clone().into() }).collect();
    let cluster = Cluster::new(events).expect("a cluster");
    let scores = HaplotypeScores::new(cluster, reads.into()).expect("reads the haplotypes explain");
    let alleles = scored
        .iter()
        .map(|allele| IndelAlleleCounts { allele: allele.clone(), total: 6 })
        .collect();
    let counts = IndelCounts { alleles };
    let scored = Some(Scored::new(Arc::new(scores), column.pos));
    column.indel_data = Some(Box::new(IndelData { counts, scored }));
    column
}
