//! Hand-built columns for the phasing tests.

use super::{Cost, FragmentId, GtAllele, PhaseCall, PhaseSet};
use crate::{
    call::{
        pair_hmm::{self, Allele, Cluster, Event, HaplotypeScores, ReadScores},
        pileup::indels::{IndelAllele, IndelAlleleCounts, IndelCounts, IndelData, Scored},
        variant_calling::{EstimatedGenotype, GenotypeTag},
    },
    metrics::{Alt, AltCall, PileupMetrics, PositionMetrics},
    sequence::{ChunkRegion, Region},
    utils::default,
};
use seqair_types::Pos1;
use seqair_types::{Base, Probability, RootMeanSquare};
use std::{num::NonZeroU8, sync::Arc};

pub(crate) use pair_hmm::fixtures::{ALT1_READ, ALT2_READ, NO, REF_READ};

pub(crate) const ALT_1: NonZeroU8 = NonZeroU8::MIN;
pub(crate) const ALT_2: NonZeroU8 = NonZeroU8::new(2).expect("2 > 0");

/// A `0/1` het's phase call in the shorthand the tests use: `set` as the VCF
/// prints it (`1001` for a block anchored on the 0-based column 1000), `first`
/// is `0` for the reference on haplotype 1, `1` for the alt, `quality` its
/// `PQ`.
pub(crate) fn call(set: u32, first: u8, quality: u8) -> PhaseCall {
    PhaseCall {
        set: PhaseSet(Pos1::new(set).expect("a 1-based phase set")),
        first: allele(first),
        quality: Cost(quality),
    }
}

pub(crate) fn mapq(mapq: u8) -> RootMeanSquare {
    [f64::from(mapq)].into_iter().collect()
}

pub(crate) fn allele(index: u8) -> GtAllele {
    NonZeroU8::new(index).map_or(GtAllele::Ref, GtAllele::Alt)
}

/// A column with whatever genotype and alts a test needs, scored by the
/// pair-HMM over `reads` against the reference and each alt's haplotype in
/// that record order (so `scored(f, [ref, alt1, alt2])` rows), and defaults
/// everywhere else. Without reads nothing was scored.
pub(crate) fn column(
    pos: u32,
    reference: Base,
    alts: &[(Base, AltCall)],
    genotype: Option<GenotypeTag>,
    reads: &[ReadScores],
) -> PileupMetrics {
    let snv_scores = (!reads.is_empty() && !alts.is_empty()).then(|| {
        let events =
            alts.iter().map(|(base, _)| Event { anchor: pos, allele: Allele::Base(*base) });
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
        snv_scores,
        phase: None,
    };
    // Mapped well enough that no site is guarded unless a test says so.
    metrics.ref_metrics.mapq = mapq(60);
    for alt in &mut metrics.alts {
        alt.metrics.mapq = mapq(60);
    }
    metrics.pos_metrics.extended.genotype = genotype.map(|genotype| EstimatedGenotype {
        genotype,
        likelihood: Probability::ONE,
        confidence: Probability::ONE,
    });
    metrics
}

/// The ordinary case: a `0/1` site whose single alt was called. Unless it has
/// no reads, it also gets [`backing`], so it passes the votes-per-allele
/// floor whatever reads the test gives it.
pub(crate) fn het(pos: u32, reference: Base, alt: Base, reads: &[ReadScores]) -> PileupMetrics {
    let mut reads = reads.to_vec();
    if !reads.is_empty() {
        reads.extend(backing(pos, [0, 1]));
    }
    column(pos, reference, &[(alt, AltCall::RealVariant)], Some(GenotypeTag::RefHet(ALT_1)), &reads)
}

/// Two confident votes for each of `haplotypes` at the column at `pos`, from
/// fragments that reach no other column and so link nothing.
pub(crate) fn backing(pos: u32, haplotypes: [usize; 2]) -> Vec<ReadScores> {
    let id = |k: u64| 1_000_000 + 4 * u64::from(pos) + k;
    haplotypes
        .into_iter()
        .enumerate()
        .flat_map(|(i, hap)| (0..2u64).map(move |copy| (id(2 * i as u64 + copy), hap)))
        .map(|(fragment, hap)| scored(fragment, explained_by(hap)))
        .collect()
}

/// A row the haplotype at `hap` explains and the other two do not.
fn explained_by(hap: usize) -> [f32; 3] {
    let mut log10 = [-4.0; 3];
    if let Some(likeliest) = log10.get_mut(hap) {
        *likeliest = -1.0;
    }
    log10
}

/// A read the pair-HMM scored: `log10` is its likelihood under the reference
/// and each scored allele's haplotype.
pub(crate) fn scored(fragment: u64, log10: [f32; 3]) -> ReadScores {
    ReadScores::new(FragmentId::new(fragment), 60, pair_hmm::row(log10))
}

/// Reads that make the pair-HMM genotype a column as a confident
/// heterozygote of `haplotypes`, and pass the indel site gate: six votes each.
/// Their fragments, from 1000 up, reach no other site, so they link nothing.
pub(crate) fn confirming(haplotypes: [usize; 2]) -> Vec<ReadScores> {
    haplotypes
        .into_iter()
        .flat_map(|hap| (0..6).map(move |i| scored(1000 + 10 * hap as u64 + i, explained_by(hap))))
        .collect()
}

pub(crate) fn insertion(bases: &[Base]) -> IndelAllele {
    IndelAllele::Insertion(bases.iter().copied().collect())
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
