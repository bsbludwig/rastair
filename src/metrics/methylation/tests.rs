use super::*;
use crate::{
    call::{
        tests::pileup::Pileup,
        variant_calling::{EstimatedGenotype, GenotypeSource, GenotypeTag},
    },
    pileups,
    sequence::Segment,
};
use seqair_types::{Base::*, Phred};

mod denovo_a_to_g;
mod denovo_other_to_c;
mod denovo_other_to_g;
mod denovo_t_to_c;
mod non_methylation;
mod original_cpg_ref_c;
mod original_cpg_ref_g;
mod proptests;

fn to_metrics(
    pileup: &Pileup,
    segment: &Segment,
    genotype: Option<EstimatedGenotype>,
) -> PileupMetrics {
    let mut metrics = pileup.metrics(segment).unwrap();
    if let Some(gt) = genotype {
        metrics.pos_metrics.genotype = Some(gt);
    }
    for alt in &mut metrics.alts {
        alt.verdict = crate::metrics::AltVerdict::Called(seqair_types::Probability::ONE);
    }
    metrics
}

#[track_caller]
fn assert_beta(methylated: Methylated, origin: CpgOrigin, expected_beta: f64) {
    let cpg = methylated
        .get(origin)
        .unwrap_or_else(|| panic!("expected {origin:?} CpG in {methylated:?}"));
    assert!(
        (*cpg.beta - expected_beta).abs() < 0.001,
        "Expected beta {expected_beta}, got {} for {origin:?}",
        *cpg.beta,
    );
}

#[track_caller]
fn assert_original(result: Methylated, expected_beta: f64) {
    assert_beta(result, CpgOrigin::Original, expected_beta);
}

#[track_caller]
fn assert_denovo(result: Methylated, expected_beta: f64) {
    assert_beta(result, CpgOrigin::DeNovo, expected_beta);
}

#[track_caller]
fn assert_none(result: Methylated) {
    assert!(result.is_empty(), "Expected no evidence, got {result:?}");
}

fn ct_genotype() -> EstimatedGenotype {
    EstimatedGenotype {
        genotype: GenotypeTag::CT,
        quality: Phred::from_phred(20),
        source: GenotypeSource::Binomial,
    }
}
