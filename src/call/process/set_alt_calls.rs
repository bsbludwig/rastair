use super::calc_ml::worth_scoring;
use crate::metrics::{AltVerdict, MetricsForAlt, PileupMetrics, Rejection};
use seqair_types::{Probability, SmallVec};
use std::num::NonZeroU32;
use tracing::instrument;

/// Call or reject each alt, and store the engine's genotype of the called
/// ones: the forests decide which are present, but an alt the engine's
/// genotype holds no copy of is not called, as an indel allele the pair-HMM
/// does not call is not (`PileupMetrics::call_engine_genotype`). An alt the
/// forests do not pass is rejected for its score, unless it was not worth
/// scoring or is a methylation candidate more than `m_max_coverage` reads
/// show.
///
/// Decided from the score alone, so a test can re-run the calls after
/// changing one.
#[instrument(level = "debug", skip_all)]
pub fn set_alt_calls(
    current: &mut PileupMetrics,
    ml_threshold: Probability,
    m_max_coverage: NonZeroU32,
) {
    let too_deep: SmallVec<bool, 4> = current
        .alts
        .iter()
        .map(|alt| {
            alt.metrics.depth > m_max_coverage.get()
                && MetricsForAlt { metrics: current, alt: &alt.metrics }.is_conversion()
        })
        .collect();
    for (alt, too_deep) in current.alts.iter_mut().zip(too_deep) {
        alt.verdict = if !worth_scoring(&current.pos_metrics, current.reference_base, alt.base) {
            AltVerdict::Rejected(Rejection::PreMl)
        } else if too_deep {
            AltVerdict::Rejected(Rejection::TooDeep(alt.verdict.ml()))
        } else {
            match alt.verdict.ml() {
                Some(ml) if ml >= ml_threshold => AltVerdict::Called(ml),
                Some(ml) => AltVerdict::Rejected(Rejection::LowScore(ml)),
                None => AltVerdict::Rejected(Rejection::NoScore),
            }
        };
    }
    current.pos_metrics.genotype = current.call_engine_genotype();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::tests::fixtures::{column, mapq};
    use seqair_types::Base::*;

    /// A column the forests would score: deep enough and mapped well.
    fn scorable(mut site: PileupMetrics) -> PileupMetrics {
        site.pos_metrics.depth = 10;
        site.pos_metrics.mapq = mapq(60);
        site
    }

    /// An alt the model passes but the engine never scored (here: no reads to
    /// score it on) has no genotype to hold it, so it is not called.
    #[test]
    fn an_alt_the_engine_did_not_score_is_not_called() {
        let mut site = scorable(column(1500, A, &[(T, true)], None, &[]));
        set_alt_calls(&mut site, Probability::new_panicky(0.4), NonZeroU32::MAX);
        let t = site.alts.first().expect("the T alt");
        assert_eq!(t.verdict, AltVerdict::Rejected(Rejection::NoGenotypeCopy(Probability::ONE)));
    }

    /// `low_ml_score` on an alt scored below the threshold, not on one too
    /// thin to score (`pre_ml`).
    #[test]
    fn an_unscored_alt_carries_pre_ml_alone() {
        let mut deep = scorable(column(1500, A, &[(T, true)], None, &[]));
        deep.find_alt_mut(T).expect("a T alt").verdict = AltVerdict::Scored(Probability::ZERO);
        let mut thin = column(1501, A, &[(T, false)], None, &[]);
        thin.pos_metrics.depth = 1;
        let threshold = Probability::new_panicky(0.4);
        set_alt_calls(&mut deep, threshold, NonZeroU32::MAX);
        set_alt_calls(&mut thin, threshold, NonZeroU32::MAX);
        let verdict = |site: &PileupMetrics| site.find_alt(T).expect("the alt").verdict;
        assert_eq!(verdict(&deep), AltVerdict::Rejected(Rejection::LowScore(Probability::ZERO)));
        assert_eq!(verdict(&thin), AltVerdict::Rejected(Rejection::PreMl));
    }
}
