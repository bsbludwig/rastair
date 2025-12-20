use crate::{
    call::methylation::ThresholdParams,
    metrics::{DenovoAdjecent, PileupMetrics},
    utils::{Base::*, IntoF64, logging::ThisIsABug},
    vcf::{InCpG, Methylated},
};
use color_eyre::{
    Result,
    eyre::{Context, ContextCompat},
};
use rastair_types::Probability;
use tracing::{Level, debug, instrument, trace, warn};

#[instrument(
    level="debug",
    skip_all,
    fields(contig = %current.contig(), pos = current.pos()),
    name = "methylation_call"
)]
pub fn call(config: &ThresholdParams, current: &PileupMetrics) -> Result<Option<Methylated>> {
    let res = call_methylation(config, current).wrap_err("Failed to call methylation")?;
    match res {
        Methylated::Unknown => Ok(None),
        _ => Ok(Some(res)),
    }
}

fn call_methylation(config: &ThresholdParams, p: &PileupMetrics) -> Result<Methylated> {
    let cpg = p.pos_metrics.cpg;
    let denovo_adj = p.pos_metrics.denovo_adj;
    let ref_base = p.ref_base();
    let sequence_context = &p.pileup.context;
    let ref_before = sequence_context.before_1;
    let ref_after = sequence_context.after_1;

    // Check for original CpG
    let original_beta = if cpg == InCpG::C || denovo_adj == DenovoAdjecent::ThisIsTheMatchingC {
        Some(ref_c(config, p)?)
    } else if cpg == InCpG::G || denovo_adj == DenovoAdjecent::ThisIsTheMatchingG {
        Some(ref_g(config, p)?)
    } else {
        None
    };

    // Check for de-novo CpG
    let denovo_beta = if p.alt(C).is_some() && ref_after == G {
        // creating new CpG
        if ref_base == T { Some(ref_t_to_c(config, p)?) } else { Some(ref_not_t_to_c(config, p)?) }
    } else if p.alt(G).is_some() && ref_before == C {
        // creating new CpG
        if ref_base == A { Some(ref_a_to_g(config, p)?) } else { Some(ref_not_a_to_g(config, p)?) }
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

fn ref_t_to_c(config: &ThresholdParams, record: &PileupMetrics) -> Result<Methylated> {
    assert_eq!(record.ref_base(), T);
    let t = &record.ref_metrics;
    let c = record.alt(C).wrap_err("Expected alt C at T->C denovo CpG site").this_is_a_bug()?;

    // T > C case: need to use strand to distinguish mod from unmod
    let c_counts = c.strand_count;
    let t_counts = t.strand_count;

    // If there's 2+ reads evidence for T on OB, assume het SNP and adjust beta
    // Note that T is the _ref_ here
    // TODO: some more sophisticated SNP calling here, taking into account baseq, mapq etc
    if t_counts.ob >= config.m_min_denovo_depth {
        // mod (reads showing T) are the ref here
        // divide by 2 assuming diploid genome
        let mod_count = t_counts.ot.f() / 2.;
        let total = c_counts.ot.f() + mod_count;
        if total > 0. {
            Ok(Methylated::DeNovoCpG { beta: Probability::new(mod_count / total).this_is_a_bug()? })
        } else {
            Ok(Methylated::NoEvidence)
        }
    } else {
        let mod_count = t_counts.ot;
        let total = c_counts.ot + t_counts.ot;
        if total > 0 {
            Ok(Methylated::DeNovoCpG {
                beta: Probability::new(mod_count.f() / total.f()).this_is_a_bug()?,
            })
        } else {
            Ok(Methylated::NoEvidence)
        }
    }
}

fn ref_not_t_to_c(_config: &ThresholdParams, record: &PileupMetrics) -> Result<Methylated> {
    assert_ne!(record.ref_base(), T);

    let t = record.alt(T);
    let c = record.alt(C).wrap_err("Expected alt C at non-T->C denovo CpG site").this_is_a_bug()?;

    // Get counts for `alt == T`, default to 0
    let t_counts = t.map(|a| a.strand_count).unwrap_or_default();

    // Ref is not T: count alt == T and alt == C separately
    let mod_count = t_counts.ot;
    let unmod = c.strand_count;
    let unmod_count = unmod.ot;

    // Check if there's evidence for T on the OB, which would be very
    // weird, ie a multi-allelic site (X->C _and_ X->T ?!)
    if t_counts.ob > 0 {
        debug!(?t_counts, "Evidence for multi-allelic SNP at het D/C site");
    }

    let total = mod_count + unmod_count;
    if total > 0 {
        Ok(Methylated::DeNovoCpG {
            beta: Probability::new(mod_count.f() / total.f()).this_is_a_bug()?,
        })
    } else {
        Ok(Methylated::NoEvidence)
    }
}

fn ref_a_to_g(config: &ThresholdParams, record: &PileupMetrics) -> Result<Methylated> {
    assert_eq!(record.ref_base(), A);
    let a = &record.ref_metrics;
    let g = record.alt(G).wrap_err("Expected alt G at A->G denovo CpG site").this_is_a_bug()?;

    // A > G case: similar logic but for OB strand
    let g_counts = g.strand_count;
    let a_counts = a.strand_count;

    // If there's 2+ reads evidence for A on OT, assume het SNP and adjust beta
    if a_counts.ot >= config.m_min_denovo_depth {
        // divide by 2 assuming diploid genome
        let mod_count = a_counts.ob.f() / 2.;
        let total = g_counts.ob.f() + mod_count;
        if total > 0. {
            Ok(Methylated::DeNovoCpG { beta: Probability::new(mod_count / total).this_is_a_bug()? })
        } else {
            Ok(Methylated::NoEvidence)
        }
    } else {
        let mod_count = a_counts.ob;
        let total = g_counts.ob + mod_count;
        if total > 0 {
            Ok(Methylated::DeNovoCpG {
                beta: Probability::new(mod_count.f() / total.f()).this_is_a_bug()?,
            })
        } else {
            Ok(Methylated::NoEvidence)
        }
    }
}

fn ref_not_a_to_g(_config: &ThresholdParams, record: &PileupMetrics) -> Result<Methylated> {
    assert_ne!(record.ref_base(), A);

    let a = record.alt(A);
    let g = record.alt(G).wrap_err("Expected alt G at non-A->G denovo CpG site").this_is_a_bug()?;

    let a_counts = a.map(|a| a.strand_count).unwrap_or_default();

    // Ref is not A: count alt == A and alt == G separately
    let mod_count = a_counts.ob;
    let unmod = g.strand_count;
    let unmod_count = unmod.ob;

    if a_counts.ot > 0 {
        debug!(?a_counts, "Evidence for multi-allelic SNP at het H/G site");
    }

    let total = mod_count + unmod_count;
    if total > 0 {
        Ok(Methylated::DeNovoCpG {
            beta: Probability::new(mod_count.f() / total.f()).this_is_a_bug()?,
        })
    } else {
        Ok(Methylated::NoEvidence)
    }
}

fn ref_c(_config: &ThresholdParams, record: &PileupMetrics) -> Result<Methylated> {
    assert_eq!(record.ref_base(), C);
    let c = &record.ref_metrics;

    // Check for non-T alternatives (possible C->N SNP)
    if tracing::enabled!(Level::TRACE)
        && *record.pos_metrics.denovo_adj
        && record.alts().iter().any(|b| *b != T)
    {
        trace!(
            pos = %record.pos(),
            "Possible C->N SNP next to a de-novo G"
        );
    }

    if let Some(t) = record.alt(T) {
        let t_counts = t.strand_count;
        let c_counts = c.strand_count;

        let mut mod_count = t_counts.ot.f();
        let unmod_count = c_counts.ot.f();

        if let Some(gt) = record.pos_metrics.genotype
            && gt.genotype.is_heterozygous()
        {
            // divide by 2 assuming diploid genome
            mod_count /= 2.;
        }

        let total = mod_count + unmod_count;
        if total > 0. {
            Ok(Methylated::OriginalCpG {
                beta: Probability::new(mod_count / total).this_is_a_bug()?,
            })
        } else {
            Ok(Methylated::NoEvidence)
        }
    } else {
        Ok(Methylated::NoEvidence)
    }
}

fn ref_g(_config: &ThresholdParams, record: &PileupMetrics) -> Result<Methylated> {
    assert_eq!(record.ref_base(), G);
    let g = &record.ref_metrics;

    // Check for non-A alternatives (possible G->N SNP)
    if tracing::enabled!(Level::TRACE)
        && *record.pos_metrics.denovo_adj
        && record.alts().iter().any(|b| *b != A)
    {
        trace!(
            pos = %record.pos(),
            "Possible G->N SNP next to a de-novo C"
        );
    }

    if let Some(a) = record.alt(A) {
        let a_counts = a.strand_count;
        let g_counts = g.strand_count;

        let mut mod_count = a_counts.ob.f();
        let unmod_count = g_counts.ob.f();

        if let Some(gt) = record.pos_metrics.genotype
            && gt.genotype.is_heterozygous()
        {
            // divide by 2 assuming diploid genome
            mod_count /= 2.;
        }

        let total = mod_count + unmod_count;
        if total > 0. {
            Ok(Methylated::OriginalCpG {
                beta: Probability::new(mod_count / total).this_is_a_bug()?,
            })
        } else {
            Ok(Methylated::NoEvidence)
        }
    } else {
        Ok(Methylated::NoEvidence)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        call::{
            pileup::Pileup,
            variant_calling::{EstimatedGenotype, GenotypeTag},
        },
        pileups,
        sequence::Segment,
    };
    use rastair_types::Probability;

    fn to_metrics(
        pileup: &Pileup,
        _segment: &Segment,
        genotype: Option<EstimatedGenotype>,
    ) -> PileupMetrics {
        let mut metrics = PileupMetrics::new(pileup.clone()).unwrap();
        if let Some(gt) = genotype {
            metrics.pos_metrics.extended.genotype = Some(gt);
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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

            assert_original_cpg(result, 0.75);
        }

        #[test]
        fn no_evidence() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [C G] OB,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

            assert_original_cpg(result, 0.75);
        }

        #[test]
        fn no_evidence() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [C G] OT,
            );
            let metrics = to_metrics(&ps[1], &seg, None);
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

            assert_no_evidence(result);
        }
    }

    mod denovo_t_to_c {
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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

            assert_denovo_cpg(result, 0.667);
        }

        #[test]
        fn no_evidence() {
            let (seg, ps) = pileups!(
                [T G] Ref,
                [T G] OB,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

            assert_none(result);
        }
    }

    mod denovo_a_to_g {
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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let metrics = to_metrics(&ps[1], &seg, None);
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

            assert_denovo_cpg(result, 0.667);
        }

        #[test]
        fn no_evidence() {
            let (seg, ps) = pileups!(
                [C A] Ref,
                [C A] OT,
            );
            let metrics = to_metrics(&ps[1], &seg, None);
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

            assert_denovo_cpg(result, 0.3);
        }

        #[test]
        fn no_evidence() {
            let (seg, ps) = pileups!(
                [A G] Ref,
                [A G] OB,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

            assert_denovo_cpg(result, 0.2);
        }

        #[test]
        fn no_evidence() {
            let (seg, ps) = pileups!(
                [C T] Ref,
                [C T] OT,
            );
            let metrics = to_metrics(&ps[1], &seg, None);
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

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
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

            assert_none(result);
        }

        #[test]
        fn wrong_context() {
            let (seg, ps) = pileups!(
                [C A] Ref,
                [C A] OT,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&ThresholdParams::default(), &metrics).unwrap();

            assert_none(result);
        }
    }
}
