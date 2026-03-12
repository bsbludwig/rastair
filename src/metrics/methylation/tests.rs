use super::*;
use crate::{
    call::{
        pileup::Pileup,
        variant_calling::{EstimatedGenotype, GenotypeTag},
    },
    metrics::AltCall,
    pileups,
    sequence::Segment,
};
use seqair_types::Probability;

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
    _segment: &Segment,
    genotype: Option<EstimatedGenotype>,
) -> PileupMetrics {
    let mut metrics = PileupMetrics::new(pileup.clone()).unwrap();
    if let Some(gt) = genotype {
        metrics.pos_metrics.extended.genotype = Some(gt);
    }
    for alt in &mut metrics.alts {
        alt.call = AltCall::RealVariant;
    }
    metrics
}

#[track_caller]
fn assert_beta(result: Option<Methylated>, origin: CpgOrigin, expected_beta: f64) {
    let methylated = result.expect("expected Some(Methylated)");
    let cpg = methylated
        .0
        .iter()
        .find(|b| b.origin == origin)
        .unwrap_or_else(|| panic!("expected {origin:?} CpG in {methylated:?}"));
    assert!(
        (*cpg.beta - expected_beta).abs() < 0.001,
        "Expected beta {expected_beta}, got {} for {origin:?}",
        *cpg.beta,
    );
}

#[track_caller]
fn assert_original(result: Option<Methylated>, expected_beta: f64) {
    assert_beta(result, CpgOrigin::Original, expected_beta);
}

#[track_caller]
fn assert_denovo(result: Option<Methylated>, expected_beta: f64) {
    assert_beta(result, CpgOrigin::DeNovo, expected_beta);
}

#[track_caller]
fn assert_none(result: Option<Methylated>) {
    assert!(
        result.as_ref().is_none_or(|m| !m.has_evidence()),
        "Expected None or no evidence, got {:?}",
        result,
    );
}

fn ct_genotype() -> EstimatedGenotype {
    EstimatedGenotype {
        genotype: GenotypeTag::CT,
        likelihood: Probability::new(0.99).unwrap(),
        confidence: Probability::new(0.99).unwrap(),
    }
}
