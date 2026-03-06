use crate::{
    call::variant_calling::GenotypeTag,
    metrics::{AltCall, DenovoAdjecent, PileupMetrics},
    utils::{Base::*, IntoF64, logging::ThisIsABug},
    vcf::{InCpG, Methylated},
};
use color_eyre::{Result, eyre::Context};
use seqair_types::{Base, Probability, Strand};
use tracing::instrument;

#[instrument(
    level="debug",
    skip_all,
    fields(contig = %current.contig(), pos = current.pos()),
    name = "methylation_call"
)]
pub fn call(current: &PileupMetrics) -> Result<Option<Methylated>> {
    let res = call_methylation(current).wrap_err("Failed to call methylation")?;
    match res {
        Methylated::Unknown => Ok(None),
        _ => Ok(Some(res)),
    }
}

fn call_methylation(p: &PileupMetrics) -> Result<Methylated> {
    let cpg = p.pos_metrics.cpg;
    let denovo_adj = p.pos_metrics.denovo_adj;
    let ref_base = p.ref_base();
    let sequence_context = &p.pileup.context;
    let ref_before = sequence_context.before_1;
    let ref_after = sequence_context.after_1;

    // Check for original CpG
    let original_beta = if cpg == InCpG::C || denovo_adj == DenovoAdjecent::ThisIsTheMatchingC {
        Some(ref_c(p)?)
    } else if cpg == InCpG::G || denovo_adj == DenovoAdjecent::ThisIsTheMatchingG {
        Some(ref_g(p)?)
    } else {
        None
    };

    // Check for de-novo CpG (only if the denovo alt is a real variant)
    let has_denovo_c = p.alts.iter().any(|a| a.base == C && a.call == AltCall::RealVariant);
    let has_denovo_g = p.alts.iter().any(|a| a.base == G && a.call == AltCall::RealVariant);
    let denovo_beta = if has_denovo_c && ref_after == G {
        // creating new CpG
        if ref_base == T { Some(ref_t_to_c(p)?) } else { Some(ref_not_t_to_c(p)?) }
    } else if has_denovo_g && ref_before == C {
        // creating new CpG
        if ref_base == A { Some(ref_a_to_g(p)?) } else { Some(ref_not_a_to_g(p)?) }
    } else {
        None
    };

    // Combine results
    match (original_beta, denovo_beta) {
        // Both original and de-novo CpG present
        (
            Some(Methylated::OriginalCpG { beta: orig }),
            Some(Methylated::DeNovoCpG { beta: denovo }),
        ) => Ok(Methylated::Both { original_beta: orig, denovo_beta: denovo }),

        // Only original CpG
        (Some(Methylated::OriginalCpG { beta: orig }), Some(Methylated::NoEvidence))
        | (Some(Methylated::OriginalCpG { beta: orig }), None) => {
            Ok(Methylated::OriginalCpG { beta: orig })
        }

        // Only de-novo CpG
        (Some(Methylated::NoEvidence), Some(Methylated::DeNovoCpG { beta: denovo }))
        | (None, Some(Methylated::DeNovoCpG { beta: denovo })) => {
            Ok(Methylated::DeNovoCpG { beta: denovo })
        }

        // No evidence from both
        (Some(Methylated::NoEvidence), Some(Methylated::NoEvidence))
        | (Some(Methylated::NoEvidence), None)
        | (None, Some(Methylated::NoEvidence)) => Ok(Methylated::NoEvidence),

        // Neither present
        (None, None) => Ok(Methylated::Unknown),

        // Unknown cases
        (Some(Methylated::Unknown), _) | (_, Some(Methylated::Unknown)) => Ok(Methylated::Unknown),

        // Shouldn't happen cases - functions returning wrong variant types
        (Some(Methylated::Both { .. }), _) | (_, Some(Methylated::Both { .. })) => {
            Err(color_eyre::eyre::eyre!("Unexpected Both variant in intermediate methylation call"))
        }
        (Some(Methylated::DeNovoCpG { .. }), Some(Methylated::OriginalCpG { .. }))
        | (Some(Methylated::DeNovoCpG { .. }), Some(Methylated::DeNovoCpG { .. }))
        | (Some(Methylated::DeNovoCpG { .. }), Some(Methylated::NoEvidence))
        | (Some(Methylated::DeNovoCpG { .. }), None)
        | (Some(Methylated::OriginalCpG { .. }), Some(Methylated::OriginalCpG { .. }))
        | (None, Some(Methylated::OriginalCpG { .. }))
        | (Some(Methylated::NoEvidence), Some(Methylated::OriginalCpG { .. })) => {
            Err(color_eyre::eyre::eyre!("Unexpected variant combination in methylation calling"))
        }
    }
}

fn ref_t_to_c(record: &PileupMetrics) -> Result<Methylated> {
    let mod_count = record.after_counts.get(Strand::OT, T, G).f();
    let unmod_count = record.after_counts.get(Strand::OT, C, G).f();
    let total = mod_count + unmod_count;
    if total == 0. {
        return Ok(Methylated::NoEvidence);
    }

    let beta = match record.pos_metrics.genotype {
        Some(gt) if gt.genotype.is_heterozygous() => calculate_het_snp_beta(mod_count, unmod_count),
        _ => mod_count / total,
    };

    Ok(Methylated::DeNovoCpG { beta: Probability::new(beta).this_is_a_bug()? })
}

fn ref_not_t_to_c(record: &PileupMetrics) -> Result<Methylated> {
    let mod_count = record.after_counts.get(Strand::OT, T, G).f();
    let unmod_count = record.after_counts.get(Strand::OT, C, G).f();
    let total = mod_count + unmod_count;
    if total == 0. {
        return Ok(Methylated::NoEvidence);
    }

    Ok(Methylated::DeNovoCpG { beta: Probability::new(mod_count / total).this_is_a_bug()? })
}

fn ref_a_to_g(record: &PileupMetrics) -> Result<Methylated> {
    let mod_count = record.before_counts.get(Strand::OB, A, C).f();
    let unmod_count = record.before_counts.get(Strand::OB, G, C).f();
    let total = mod_count + unmod_count;
    if total == 0. {
        return Ok(Methylated::NoEvidence);
    }

    let beta = match record.pos_metrics.genotype {
        Some(gt) if gt.genotype.is_heterozygous() => calculate_het_snp_beta(mod_count, unmod_count),
        _ => mod_count / total,
    };

    Ok(Methylated::DeNovoCpG { beta: Probability::new(beta).this_is_a_bug()? })
}

fn ref_not_a_to_g(record: &PileupMetrics) -> Result<Methylated> {
    let mod_count = record.before_counts.get(Strand::OB, A, C).f();
    let unmod_count = record.before_counts.get(Strand::OB, G, C).f();
    let total = mod_count + unmod_count;
    if total == 0. {
        return Ok(Methylated::NoEvidence);
    }

    Ok(Methylated::DeNovoCpG { beta: Probability::new(mod_count / total).this_is_a_bug()? })
}

/// Calculate beta value for heterozygous C>T or G>A variants.
///
/// For heterozygous SNPs, assume 50% of the modified reads are expected to come
/// from the SNP allele itself, not from methylation. So we calculate beta by
/// only counting the "excess" modified reads.
fn calculate_het_snp_beta(mod_count: f64, unmod_count: f64) -> f64 {
    let total = mod_count + unmod_count;
    if total == 0. {
        return 0.;
    }

    let excess_mod = (mod_count - total / 2.).max(0.0); // mod_count * (0.5 - unmod_fraction);
    excess_mod / (unmod_count + excess_mod)
}

fn het_alt_is_base(record: &PileupMetrics, base: Base) -> bool {
    let Some(gt) = record.pos_metrics.genotype else {
        return false;
    };

    match gt.genotype {
        GenotypeTag::RefHet(idx) => {
            let i = (idx.get() as usize).saturating_sub(1);
            record.alts.get(i).map(|a| a.base) == Some(base)
        }
        GenotypeTag::AltHet(a, b) => {
            let ai = (a.get() as usize).saturating_sub(1);
            let bi = (b.get() as usize).saturating_sub(1);
            record.alts.get(ai).map(|x| x.base) == Some(base)
                || record.alts.get(bi).map(|x| x.base) == Some(base)
        }
        _ => false,
    }
}

fn ref_c(record: &PileupMetrics) -> Result<Methylated> {
    let mod_count = record.after_counts.get(Strand::OT, T, G).f();
    let unmod_count = record.after_counts.get(Strand::OT, C, G).f();
    let total = mod_count + unmod_count;
    if total == 0. {
        return Ok(Methylated::NoEvidence);
    }

    let beta_value = match record.pos_metrics.genotype {
        Some(gt) if gt.genotype.is_heterozygous() && het_alt_is_base(record, T) => {
            calculate_het_snp_beta(mod_count, unmod_count)
        }
        Some(gt) if gt.genotype.is_homozygous() && !gt.genotype.is_hom_ref() => 0.0,
        _ => mod_count / total,
    };

    Ok(Methylated::OriginalCpG { beta: Probability::new(beta_value).this_is_a_bug()? })
}

fn ref_g(record: &PileupMetrics) -> Result<Methylated> {
    let mod_count = record.before_counts.get(Strand::OB, A, C).f();
    let unmod_count = record.before_counts.get(Strand::OB, G, C).f();
    let total = mod_count + unmod_count;
    if total == 0. {
        return Ok(Methylated::NoEvidence);
    }

    let beta_value = match record.pos_metrics.genotype {
        Some(gt) if gt.genotype.is_heterozygous() && het_alt_is_base(record, A) => {
            calculate_het_snp_beta(mod_count, unmod_count)
        }
        Some(gt) if gt.genotype.is_homozygous() && !gt.genotype.is_hom_ref() => 0.0,
        _ => mod_count / total,
    };

    Ok(Methylated::OriginalCpG { beta: Probability::new(beta_value).this_is_a_bug()? })
}

#[cfg(test)]
mod tests {
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

    fn to_metrics(
        pileup: &Pileup,
        _segment: &Segment,
        genotype: Option<EstimatedGenotype>,
    ) -> PileupMetrics {
        let mut metrics = PileupMetrics::new(pileup.clone()).unwrap();
        if let Some(gt) = genotype {
            metrics.pos_metrics.extended.genotype = Some(gt);
        }
        // Methylation calling now respects alt calls for de-novo detection, so
        // mark all observed alts as real variants in these unit tests.
        for alt in &mut metrics.alts {
            alt.call = AltCall::RealVariant;
        }
        metrics
    }

    #[track_caller]
    fn assert_original_cpg(result: Option<Methylated>, expected_beta: f64) {
        match result {
            Some(Methylated::OriginalCpG { beta }) => {
                assert!(
                    (*beta - expected_beta).abs() < 0.001,
                    "Expected beta {}, got {}",
                    expected_beta,
                    *beta
                );
            }
            other => panic!("Expected OriginalCpG, got {:?}", other),
        }
    }

    #[track_caller]
    fn assert_denovo_cpg(result: Option<Methylated>, expected_beta: f64) {
        match result {
            Some(Methylated::DeNovoCpG { beta }) => {
                assert!(
                    (*beta - expected_beta).abs() < 0.001,
                    "Expected beta {}, got {}",
                    expected_beta,
                    *beta
                );
            }
            other => panic!("Expected DeNovoCpG, got {:?}", other),
        }
    }

    #[track_caller]
    fn assert_no_evidence(result: Option<Methylated>) {
        match result {
            Some(Methylated::NoEvidence) => {}
            other => panic!("Expected NoEvidence, got {:?}", other),
        }
    }

    #[track_caller]
    fn assert_none(result: Option<Methylated>) {
        assert!(result.is_none(), "Expected None (Unknown), got {:?}", result);
    }

    fn ct_genotype() -> EstimatedGenotype {
        EstimatedGenotype {
            genotype: GenotypeTag::CT,
            likelihood: Probability::new(0.99).unwrap(),
            confidence: Probability::new(0.99).unwrap(),
        }
    }

    mod original_cpg_ref_c {
        use super::*;

        #[test]
        fn fully_methylated() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [T A] OT,
                [T A] OT,
                [T A] OT,
                [C G] OB,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&metrics).unwrap();

            assert_original_cpg(result, 1.0);
        }

        #[test]
        fn unmethylated() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [C G] OT,
                [C G] OT,
                [C G] OT,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&metrics).unwrap();

            assert_no_evidence(result);
        }

        #[test]
        fn partially_methylated() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [T A] OT,
                [T A] OT,
                [T A] OT,
                [C G] OT,
                [C G] OT,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&metrics).unwrap();

            assert_original_cpg(result, 0.6);
        }

        #[test]
        fn het_snp() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [T A] OT,
                [T A] OT,
                [T A] OT,
                [T A] OT,
                [T A] OT,
                [T A] OT,
                [C G] OT,
            );

            let metrics = to_metrics(&ps[0], &seg, Some(ct_genotype()));
            let result = call(&metrics).unwrap();

            assert_original_cpg(result, 0.71428571);
        }

        #[test]
        fn het_snp_fifty_fifty() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [T A] OT,
                [T A] OT,
                [T A] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
            );
            let metrics = to_metrics(&ps[0], &seg, Some(ct_genotype()));
            let result = call(&metrics).unwrap();

            // With corrected formula: 3 mod, 3 unmod (50/50 split)
            // excess_mod = 3 * (0.5 - 3/6) = 3 * 0 = 0
            // beta = 0
            assert_original_cpg(result, 0.0);
        }

        #[test]
        fn no_evidence() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [C G] OB,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&metrics).unwrap();

            assert_no_evidence(result);
        }

        #[test]
        fn no_alt_t() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [C G] OT,
                [C G] OT,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&metrics).unwrap();

            assert_no_evidence(result);
        }
    }

    mod original_cpg_ref_g {
        use super::*;

        #[test]
        fn fully_methylated() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [C A] OB,
                [C A] OB,
                [C A] OB,
                [C G] OT,
            );
            let metrics = to_metrics(&ps[1], &seg, None);
            let result = call(&metrics).unwrap();

            assert_original_cpg(result, 1.0);
        }

        #[test]
        fn unmethylated() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [C G] OB,
                [C G] OB,
                [C G] OB,
            );
            let metrics = to_metrics(&ps[1], &seg, None);
            let result = call(&metrics).unwrap();

            assert_no_evidence(result);
        }

        #[test]
        fn partially_methylated() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [C A] OB,
                [C A] OB,
                [C A] OB,
                [C A] OB,
                [C A] OB,
                [C A] OB,
                [C A] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
            );
            let metrics = to_metrics(&ps[1], &seg, None);
            let result = call(&metrics).unwrap();

            assert_original_cpg(result, 0.7);
        }

        #[test]
        fn het_snp() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [C A] OB,
                [C A] OB,
                [C A] OB,
                [C A] OB,
                [C A] OB,
                [C A] OB,
                [C G] OB,
            );
            let metrics = to_metrics(&ps[1], &seg, Some(ct_genotype()));
            let result = call(&metrics).unwrap();

            assert_original_cpg(result, 0.71428571);
        }

        #[test]
        fn het_snp_fifty_fifty() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [C A] OB,
                [C A] OB,
                [C A] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
            );
            let metrics = to_metrics(&ps[1], &seg, Some(ct_genotype()));
            let result = call(&metrics).unwrap();

            // With corrected formula: 3 mod, 3 unmod (50/50 split)
            // excess_mod = 3 * (0.5 - 3/6) = 3 * 0 = 0
            // beta = 0
            assert_original_cpg(result, 0.0);
        }

        #[test]
        fn no_evidence() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [C G] OT,
            );
            let metrics = to_metrics(&ps[1], &seg, None);
            let result = call(&metrics).unwrap();

            assert_no_evidence(result);
        }

        #[test]
        fn no_alt_a() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [C G] OB,
                [C G] OB,
            );
            let metrics = to_metrics(&ps[1], &seg, None);
            let result = call(&metrics).unwrap();

            assert_no_evidence(result);
        }

        #[test]
        fn filtered_denovo_alt_does_not_change_beta() {
            // CGG context with a real G>T variant and a filtered-out G>C denovo candidate.
            // The denovo beta must not be used when the denovo alt is not a real variant.
            let (seg, ps) = pileups!(
                [C G G] Ref,
                [C T G] OT,
                [C T G] OT,
                [C C G] OB,
            );
            let mut metrics = to_metrics(&ps[1], &seg, None);

            let denovo_alt = metrics.alts.iter_mut().find(|a| a.base == C).expect("expected C alt");
            denovo_alt.call = AltCall::ReadError;

            let result = call(&metrics).unwrap();
            assert_no_evidence(result);
        }
    }

    mod denovo_t_to_c {
        use std::num::NonZeroU8;

        use super::*;

        #[test]
        fn standard() {
            let (seg, ps) = pileups!(
                [T G] Ref,
                [T G] OT,
                [T G] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&metrics).unwrap();

            assert_denovo_cpg(result, 0.2);
        }

        #[test]
        fn het_snp_adjustment() {
            let (seg, ps) = pileups!(
                [T G] Ref,
                [T G] OT,
                [T G] OT,
                [T G] OT,
                [T G] OT,
                [C G] OT,
                [T G] OB,
                [T G] OB,
            );
            let gt = EstimatedGenotype {
                genotype: GenotypeTag::RefHet(NonZeroU8::new(1).unwrap()),
                likelihood: Probability::new(0.8).unwrap(),
                confidence: Probability::new(0.99).unwrap(),
            };

            let metrics = to_metrics(&ps[0], &seg, Some(gt));
            let result = call(&metrics).unwrap();

            assert_denovo_cpg(result, 0.6);
        }

        #[test]
        fn no_evidence() {
            let (seg, ps) = pileups!(
                [T G] Ref,
                [T G] OB,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&metrics).unwrap();

            assert_none(result);
        }
    }

    mod denovo_a_to_g {
        use std::num::NonZeroU8;

        use super::*;

        #[test]
        fn standard() {
            let (seg, ps) = pileups!(
                [C A] Ref,
                [C A] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
            );
            let metrics = to_metrics(&ps[1], &seg, None);
            let result = call(&metrics).unwrap();

            assert_denovo_cpg(result, 0.1);
        }

        #[test]
        fn het_snp_adjustment() {
            let (seg, ps) = pileups!(
                [C A] Ref,
                [C A] OB,
                [C A] OB,
                [C A] OB,
                [C A] OB,
                [C G] OB,
                [C A] OT,
                [C A] OT,
            );
            let gt = EstimatedGenotype {
                genotype: GenotypeTag::RefHet(NonZeroU8::new(1).unwrap()),
                likelihood: Probability::new(0.8).unwrap(),
                confidence: Probability::new(0.99).unwrap(),
            };
            let metrics = to_metrics(&ps[1], &seg, Some(gt));
            let result = call(&metrics).unwrap();

            assert_denovo_cpg(result, 0.6);
        }

        #[test]
        fn no_evidence() {
            let (seg, ps) = pileups!(
                [C A] Ref,
                [C A] OT,
            );
            let metrics = to_metrics(&ps[1], &seg, None);
            let result = call(&metrics).unwrap();

            assert_none(result);
        }
    }

    mod denovo_other_to_c {
        use super::*;

        #[test]
        fn standard_a_to_c() {
            let (seg, ps) = pileups!(
                [A G] Ref,
                [T G] OT,
                [T G] OT,
                [T G] OT,
                [T G] OT,
                [T G] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
            );
            let gt = EstimatedGenotype {
                genotype: GenotypeTag::HomRef,
                likelihood: Probability::new(0.8).unwrap(),
                confidence: Probability::new(0.99).unwrap(),
            };
            let metrics = to_metrics(&ps[0], &seg, Some(gt));
            let result = call(&metrics).unwrap();

            assert_denovo_cpg(result, 0.5);
        }

        #[test]
        fn multi_allelic_warning() {
            let (seg, ps) = pileups!(
                [A G] Ref,
                [T G] OT,
                [T G] OT,
                [T G] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
                [T G] OB,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&metrics).unwrap();

            assert_denovo_cpg(result, 0.3);
        }

        #[test]
        fn no_evidence() {
            let (seg, ps) = pileups!(
                [A G] Ref,
                [A G] OB,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&metrics).unwrap();

            assert_none(result);
        }
    }

    mod denovo_other_to_g {
        use super::*;

        #[test]
        fn standard_t_to_g() {
            let (seg, ps) = pileups!(
                [C T] Ref,
                [C A] OB,
                [C A] OB,
                [C A] OB,
                [C A] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
            );
            let metrics = to_metrics(&ps[1], &seg, None);
            let result = call(&metrics).unwrap();

            assert_denovo_cpg(result, 0.4);
        }

        #[test]
        fn multi_allelic_warning() {
            let (seg, ps) = pileups!(
                [C T] Ref,
                [C A] OB,
                [C A] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
                [C G] OB,
                [C A] OT,
            );
            let metrics = to_metrics(&ps[1], &seg, None);
            let result = call(&metrics).unwrap();

            assert_denovo_cpg(result, 0.2);
        }

        #[test]
        fn no_evidence() {
            let (seg, ps) = pileups!(
                [C T] Ref,
                [C T] OT,
            );
            let metrics = to_metrics(&ps[1], &seg, None);
            let result = call(&metrics).unwrap();

            assert_none(result);
        }
    }

    mod non_methylation {
        use super::*;

        #[test]
        fn non_cpg_context() {
            let (seg, ps) = pileups!(
                [A T] Ref,
                [A T] OT,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&metrics).unwrap();

            assert_none(result);
        }

        #[test]
        fn wrong_context() {
            let (seg, ps) = pileups!(
                [C A] Ref,
                [C A] OT,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&metrics).unwrap();

            assert_none(result);
        }
    }
}
