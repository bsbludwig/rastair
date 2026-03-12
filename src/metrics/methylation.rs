use crate::{
    call::variant_calling::GenotypeTag,
    metrics::{AltCall, DenovoAdjecent, PileupMetrics, ReadKey},
    utils::{Base::*, IntoF64, logging::ThisIsABug},
    vcf::{CpgBeta, CpgOrigin, InCpG, Methylated},
};
use color_eyre::{Result, eyre::Context};
use seqair_types::{Base, Probability, SmallVec, Strand};
use tracing::instrument;

#[instrument(
    level="debug",
    skip_all,
    fields(contig = %current.contig(), pos = current.pos()),
    name = "methylation_call"
)]
pub fn call(current: &PileupMetrics) -> Result<Option<Methylated>> {
    let mut betas: SmallVec<CpgBeta, 2> = SmallVec::new();

    if let Some(b) = compute_beta(current, CpgSide::C).wrap_err("C-side beta")? {
        betas.push(b);
    }
    if let Some(b) = compute_beta(current, CpgSide::G).wrap_err("G-side beta")? {
        betas.push(b);
    }

    Ok((!betas.is_empty()).then_some(Methylated(betas)))
}

/// Determine whether this position has a CpG allele on the given side,
/// and if so compute the methylation beta value.
fn compute_beta(record: &PileupMetrics, side: CpgSide) -> Result<Option<CpgBeta>> {
    let Some(origin) = cpg_origin(record, side) else { return Ok(None) };

    let (raw_mod, raw_unmod) = read_counts(record, side);
    let mod_count = raw_mod.f();
    let unmod_count = raw_unmod.f();

    if mod_count + unmod_count == 0. {
        return Ok(None);
    }

    let adjustment = genotype_adjustment(record, side, origin);
    let beta = adjusted_beta(mod_count, unmod_count, adjustment);

    Ok(Some(CpgBeta {
        origin,
        beta: Probability::new(beta).this_is_a_bug()?,
        mod_count: raw_mod,
        total_count: raw_mod + raw_unmod,
    }))
}

/// Which side of the CpG dinucleotide we are looking at.
#[derive(Debug, Clone, Copy)]
enum CpgSide {
    C,
    G,
}

impl CpgSide {
    fn strand(self) -> Strand {
        match self {
            CpgSide::C => Strand::OT,
            CpgSide::G => Strand::OB,
        }
    }

    /// The base that appears when methylated (T for C-side, A for G-side).
    fn mod_base(self) -> Base {
        match self {
            CpgSide::C => T,
            CpgSide::G => A,
        }
    }

    /// The base that appears when unmethylated (C for C-side, G for G-side).
    fn unmod_base(self) -> Base {
        match self {
            CpgSide::C => C,
            CpgSide::G => G,
        }
    }

    /// The required adjacent base (G after C, C before G).
    fn adjacent_base(self) -> Base {
        match self {
            CpgSide::C => G,
            CpgSide::G => C,
        }
    }
}

/// Determine the CpG origin for this position on the given side, or None if
/// this position doesn't have a CpG allele on that side.
fn cpg_origin(record: &PileupMetrics, side: CpgSide) -> Option<CpgOrigin> {
    let is_original = match side {
        CpgSide::C => {
            record.pos_metrics.cpg == InCpG::C
                || record.pos_metrics.denovo_adj == DenovoAdjecent::ThisIsTheMatchingC
        }
        CpgSide::G => {
            record.pos_metrics.cpg == InCpG::G
                || record.pos_metrics.denovo_adj == DenovoAdjecent::ThisIsTheMatchingG
        }
    };

    if is_original {
        return Some(CpgOrigin::Original);
    }

    let cpg_base = side.unmod_base();
    let has_denovo_alt =
        record.alts.iter().any(|a| a.base == cpg_base && a.call == AltCall::RealVariant);

    let adjacent_present = match side {
        CpgSide::C => record.pileup.context.after_1 == Some(G),
        CpgSide::G => record.pileup.context.before_1 == Some(C),
    };

    if has_denovo_alt && adjacent_present {
        return Some(CpgOrigin::DeNovo);
    }

    None
}

/// Look up the mod and unmod read counts for the given CpG side.
fn read_counts(record: &PileupMetrics, side: CpgSide) -> (u32, u32) {
    let strand = side.strand();
    let adj = side.adjacent_base();

    let (counts, mod_base, unmod_base) = match side {
        CpgSide::C => (&record.after_counts, T, C),
        CpgSide::G => (&record.before_counts, A, G),
    };

    let raw_mod = counts.get(ReadKey { strand, current: mod_base, adj });
    let raw_unmod = counts.get(ReadKey { strand, current: unmod_base, adj });
    (raw_mod, raw_unmod)
}

/// How the genotype at this position affects the beta calculation.
enum GenotypeAdjustment {
    /// Normal beta: mod / (mod + unmod).
    None,
    /// The confounding base (T for C-side, A for G-side) is present as a het
    /// allele, so some mod-base reads come from the SNP, not methylation.
    HetConfounded,
    /// Original ref base is fully replaced by a homozygous variant — the
    /// original CpG no longer exists on either chromosome.
    HomAlt,
}

fn genotype_adjustment(
    record: &PileupMetrics,
    side: CpgSide,
    origin: CpgOrigin,
) -> GenotypeAdjustment {
    let Some(gt) = record.pos_metrics.genotype else {
        return GenotypeAdjustment::None;
    };

    if origin == CpgOrigin::Original && gt.genotype.is_homozygous() && !gt.genotype.is_hom_ref() {
        return GenotypeAdjustment::HomAlt;
    }

    if gt.genotype.is_heterozygous() {
        let confounding = side.mod_base();
        let confounded = match origin {
            CpgOrigin::Original => het_alt_is_base(record, &gt.genotype, confounding),
            CpgOrigin::DeNovo => {
                record.ref_base() == confounding
                    || het_alt_is_base(record, &gt.genotype, confounding)
            }
        };
        if confounded {
            return GenotypeAdjustment::HetConfounded;
        }
    }

    GenotypeAdjustment::None
}

fn het_alt_is_base(record: &PileupMetrics, gt: &GenotypeTag, base: Base) -> bool {
    match *gt {
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

fn adjusted_beta(mod_count: f64, unmod_count: f64, adjustment: GenotypeAdjustment) -> f64 {
    match adjustment {
        GenotypeAdjustment::HomAlt => 0.0,
        GenotypeAdjustment::HetConfounded => {
            let total = mod_count + unmod_count;
            let excess_mod = (mod_count - total / 2.).max(0.0);
            excess_mod / (unmod_count + excess_mod)
        }
        GenotypeAdjustment::None => mod_count / (mod_count + unmod_count),
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

    mod original_cpg_ref_c {
        use super::*;

        #[test]
        fn fully_methylated() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [T G] OT,
                [T G] OT,
                [T G] OT,
                [C G] OB,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&metrics).unwrap();

            assert_original(result, 1.0);
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

            assert_original(result, 0.0);
        }

        #[test]
        fn partially_methylated() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [T G] OT,
                [T G] OT,
                [T G] OT,
                [C G] OT,
                [C G] OT,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&metrics).unwrap();

            assert_original(result, 0.6);
        }

        #[test]
        fn het_snp() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [T G] OT,
                [T G] OT,
                [T G] OT,
                [T G] OT,
                [T G] OT,
                [T G] OT,
                [C G] OT,
            );

            let metrics = to_metrics(&ps[0], &seg, Some(ct_genotype()));
            let result = call(&metrics).unwrap();

            assert_original(result, 0.71428571);
        }

        #[test]
        fn het_snp_fifty_fifty() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [T G] OT,
                [T G] OT,
                [T G] OT,
                [C G] OT,
                [C G] OT,
                [C G] OT,
            );
            let metrics = to_metrics(&ps[0], &seg, Some(ct_genotype()));
            let result = call(&metrics).unwrap();

            assert_original(result, 0.0);
        }

        #[test]
        fn no_evidence() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [C G] OB,
            );
            let metrics = to_metrics(&ps[0], &seg, None);
            let result = call(&metrics).unwrap();

            assert_none(result);
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

            assert_original(result, 0.0);
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

            assert_original(result, 1.0);
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

            assert_original(result, 0.0);
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

            assert_original(result, 0.7);
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

            assert_original(result, 0.71428571);
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

            assert_original(result, 0.0);
        }

        #[test]
        fn no_evidence() {
            let (seg, ps) = pileups!(
                [C G] Ref,
                [C G] OT,
            );
            let metrics = to_metrics(&ps[1], &seg, None);
            let result = call(&metrics).unwrap();

            assert_none(result);
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

            assert_original(result, 0.0);
        }

        #[test]
        fn filtered_denovo_alt_does_not_change_beta() {
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
            assert_none(result);
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

            assert_denovo(result, 0.2);
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

            assert_denovo(result, 0.6);
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

            assert_denovo(result, 0.1);
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

            assert_denovo(result, 0.6);
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

            assert_denovo(result, 0.5);
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

            assert_denovo(result, 0.3);
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

            assert_denovo(result, 0.4);
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

            assert_denovo(result, 0.2);
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
