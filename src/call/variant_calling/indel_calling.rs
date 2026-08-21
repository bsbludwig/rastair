use super::GenotypeTag;
use crate::call::pileup::indels::{IndelAllele, IndelCounts};
use better_default::Default;
use probability::prelude::{Binomial, Discrete as _};
use seqair_types::{Phred, Probability};
use std::num::NonZeroU8;
use tracing::{instrument, trace};

/// CLI parameters for indel calling.
#[derive(Debug, Clone, clap::Args, serde::Serialize, serde::Deserialize, Default)]
pub struct IndelParams {
    /// Enable experimental indel calling using hard filters
    ///
    /// Candidate indels are accepted by a fixed filter chain with no ML scoring.
    /// When disabled, Rastair calls SNPs and methylation only.
    #[arg(long, default_value_t = false)]
    pub experimental_indels: bool,

    /// Enable experimental indel calling scored by the ML model
    ///
    /// Selects the ML indel pathway instead of the hard-filter chain. Unlike
    /// `--no-ml` this leaves the SNV and CpG models alone, so the indel pathway can
    /// be chosen independently of them. Degrades to the hard-filter chain under
    /// `--no-ml` rather than erroring.
    #[arg(long, default_value_t = false, conflicts_with = "experimental_indels")]
    pub experimental_indels_ml: bool,

    /// Let the model rescue indels the binomial genotyped homozygous reference.
    ///
    /// The two pathways currently intersect: `--experimental-indels-ml` lets a low
    /// score veto an allele the hard chain would have kept, which measured costs
    /// more recall than it buys precision. This unions them instead -- the hard
    /// chain's calls stand, and an allele it dropped as hom-ref is reinstated as a
    /// heterozygote when the model rates it above the threshold.
    #[arg(long, default_value_t = false, conflicts_with = "experimental_indels_ml")]
    pub indel_ml_rescue: bool,

    /// Minimum alternate observations to call an indel
    ///
    /// This threshold is evaluated per indel allele after read-level indel
    /// filters are applied.
    #[arg(long, default_value_t = 2)]
    #[arg(help_heading = crate::utils::cli::sections::FILTER)]
    #[default(2)]
    pub min_indel_ao: u32,

    /// Minimum depth to call an indel
    ///
    /// Depth is computed at the locus after applying indel-specific read
    /// filters and depth adjustments.
    #[arg(long, default_value_t = 2)]
    #[arg(help_heading = crate::utils::cli::sections::FILTER)]
    #[default(2)]
    pub min_indel_depth: u32,

    /// Error rate for indel genotyping (higher than SNV due to alignment uncertainty)
    #[arg(long, default_value_t = 0.05)]
    #[arg(help_heading = crate::utils::cli::sections::PROCESSING)]
    #[default(0.05)]
    pub indel_error_rate: f64,

    /// Maximum number of non-TAPS mismatches allowed on a read supporting an indel.
    ///
    /// C->T mismatches on OT reads and G->A mismatches on OB reads are excluded
    /// from the count as they are expected TAPS methylation signal.
    #[arg(long, default_value_t = 5)]
    #[arg(help_heading = crate::utils::cli::sections::FILTER)]
    #[default(5)]
    pub indel_max_mismatches: u32,

    /// Ignore indels within this many bases of either read end.
    ///
    /// This stricter end-of-read filter helps suppress alignment artifacts at
    /// read starts and ends.
    #[arg(long, default_value_t = 0)]
    #[arg(help_heading = crate::utils::cli::sections::FILTER)]
    #[default(0)]
    pub indel_end_of_read_cutoff: usize,
}

/// Result of indel calling at a single position for one allele.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IndelCall {
    pub allele: IndelAllele,
    pub genotype: GenotypeTag,
    pub quality: Phred,
    /// ML score: probability this indel is a true variant.
    pub ml: Option<Probability>,
    /// Depth with noisy reads removed from both sides; see [`IndelCounts::clean_depth`].
    pub depth: u32,
    /// Reads supporting this indel allele.
    pub alt_count: u32,
    /// Supported on only one bisulfite strand. Emitted with an `indel_strand`
    /// FILTER rather than dropped, so `--all` shows why the allele did not pass.
    pub one_sided: bool,
}

/// Call indels at a position. Returns empty vec if no indels pass filters.
///
/// When `ml_enabled` is true, the binomial genotype test is used only for
/// informational genotyping — all alleles passing min AO and min depth
/// are forwarded to the ML model. When ML is off, the binomial test
/// acts as a hard gate (`hom_ref` alleles are rejected).
#[instrument(level = "trace", skip_all)]
pub fn call_indels(
    indels: &IndelCounts,
    params: &IndelParams,
    ml_enabled: bool,
    tract: u32,
    rescue: bool,
) -> Vec<IndelCall> {
    let mut calls = Vec::new();

    if indels.is_empty() {
        return calls;
    }

    let filtered_depth = indels.clean_depth();

    if filtered_depth < params.min_indel_depth {
        return calls;
    }

    for allele_counts in &indels.alleles {
        let alt_count = allele_counts.clean_total();
        if alt_count < params.min_indel_ao {
            trace!(
                allele = ?allele_counts.allele,
                alt_count,
                min = params.min_indel_ao,
                "Indel skipped: below min AO"
            );
            continue;
        }

        let genotype = binomial_genotype(
            alt_count as usize,
            filtered_depth as usize,
            params.indel_error_rate,
            tract,
        );

        if !ml_enabled {
            let Some(g) = genotype else { continue };
            // With rescue on, a hom-ref allele is kept so the model can reinstate it
            // after scoring; `rescue_hom_ref` drops the ones it does not.
            if matches!(g.tag, GenotypeTag::HomRef) && !rescue {
                trace!(
                    allele = ?allele_counts.allele,
                    alt_count,
                    depth = filtered_depth,
                    "Indel skipped: genotyped as hom ref (ML off)"
                );
                continue;
            }
            calls.push(IndelCall {
                allele: allele_counts.allele.clone(),
                genotype: g.tag,
                quality: g.quality,
                ml: None,
                depth: filtered_depth,
                alt_count,
                one_sided: !allele_counts.on_both_strands(),
            });
        } else {
            let (tag, quality) = genotype
                .map(|g| (g.tag, g.quality))
                .unwrap_or((GenotypeTag::hom_ref(), Phred::from_phred(0_u8)));
            calls.push(IndelCall {
                allele: allele_counts.allele.clone(),
                genotype: tag,
                quality,
                ml: None,
                depth: filtered_depth,
                alt_count,
                one_sided: !allele_counts.on_both_strands(),
            });
        }
    }

    resolve_compound_het(&mut calls, filtered_depth, params.indel_error_rate, tract);
    calls
}

impl IndelParams {
    pub fn enabled(&self) -> bool {
        self.experimental_indels || self.experimental_indels_ml || self.indel_ml_rescue
    }

    /// Whether indels need ML scores at all: either the model decides, or it is
    /// there to rescue what the binomial rejected.
    pub fn needs_ml_scores(&self, ml_enabled: bool) -> bool {
        (self.experimental_indels_ml || self.indel_ml_rescue) && ml_enabled
    }

    /// Whether the ML pathway is selected *and* a model is available.
    pub fn use_ml(&self, ml_enabled: bool) -> bool {
        self.experimental_indels_ml && ml_enabled
    }
}

/// Measured allele fraction of a real indel, by the tandem tract it sits in:
/// `(tract bases, het, hom-alt)`. Fitted on ~60k GIAB-judged indels across a 5base
/// and a TAPS sample; interpolated between rows, held outside them.
///
/// Reference bias keeps a hom-alt indel near 0.90 even in simple sequence, and past
/// a tract of ~11 it falls to 0.68 while het only falls to 0.37 — so a fixed
/// `1 - error_rate` hypothesis rejects almost every real hom-alt call.
const ALLELE_FRACTION_BY_TRACT: [(f64, f64, f64); 5] = [
    (6.0, 0.440, 0.900),
    (9.0, 0.439, 0.883),
    (12.5, 0.409, 0.779),
    (17.5, 0.379, 0.683),
    (25.0, 0.366, 0.683),
];

/// Measured *pooled* fraction of a compound heterozygote, by tandem tract. A `1/2`
/// locus has no reference chromosome but still presents at ~0.56 in simple
/// sequence, nothing like a hom-alt.
const COMPOUND_FRACTION_BY_TRACT: [(f64, f64); 5] =
    [(6.0, 0.560), (9.0, 0.780), (12.5, 0.759), (17.5, 0.660), (25.0, 0.619)];

fn interpolate<const N: usize>(table: [(f64, f64); N], tract: u32) -> f64 {
    let tract = f64::from(tract);
    if tract <= table[0].0 {
        return table[0].1;
    }
    for pair in table.windows(2) {
        let [low, high] = pair else { continue };
        if tract <= high.0 {
            return low.1 + (tract - low.0) / (high.0 - low.0) * (high.1 - low.1);
        }
    }
    table[N - 1].1
}

/// Expected `(het, hom_alt)` allele fraction in a tract of `tract` bases.
fn expected_allele_fractions(tract: u32) -> (f64, f64) {
    let het = ALLELE_FRACTION_BY_TRACT.map(|(t, h, _)| (t, h));
    let hom = ALLELE_FRACTION_BY_TRACT.map(|(t, _, h)| (t, h));
    (interpolate(het, tract), interpolate(hom, tract))
}

/// Reinstate hom-ref indel calls the model rates as real, and drop the rest.
///
/// Runs after ML scoring, which is the first point where both the binomial's
/// verdict and the model's score are available. A rescued allele becomes a
/// heterozygote: it is the most conservative genotype that is still a variant, and
/// the binomial having preferred hom-ref means the evidence is at the low end.
pub fn rescue_hom_ref(calls: &mut Vec<IndelCall>, ml_threshold: Option<Probability>) {
    let alt_one = NonZeroU8::new(1).expect("nonzero");
    calls.retain_mut(|call| {
        if !matches!(call.genotype, GenotypeTag::HomRef) {
            // The score is dropped on calls the binomial already accepted: this mode
            // exists so the model can *add* calls, and leaving a score attached would
            // let `low_ml_score` veto them, which is the behaviour being replaced.
            call.ml = None;
            return true;
        }
        match (call.ml, ml_threshold) {
            (Some(ml), Some(threshold)) if ml >= threshold => {
                call.genotype = GenotypeTag::ref_het(alt_one);
                true
            }
            _ => false,
        }
    });
}

/// Promote a locus whose two best alleles account for it together, and where
/// neither does so alone, to a compound heterozygote.
///
/// Weighed against the pooled count: one real allele plus noise, two real alleles
/// and no reference chromosome, or one allele carrying the locus. Each allele must
/// also stand on its own, because a compound het's pooled fraction is close to a
/// het carrying a trace second allele.
fn resolve_compound_het(calls: &mut Vec<IndelCall>, depth: u32, error_rate: f64, tract: u32) {
    if calls.len() < 2 || depth == 0 {
        return;
    }
    // A one-sided allele is emitted for visibility but is not evidence, so it must
    // not consume one of the two compound-het slots.
    calls.sort_by_key(|c| (c.one_sided, std::cmp::Reverse(c.alt_count)));
    let (Some(first), Some(second)) = (calls.first(), calls.get(1)) else { return };
    if first.one_sided || second.one_sided {
        return;
    }
    if matches!(first.genotype, GenotypeTag::HomAlt(_)) {
        return;
    }
    let stands_alone = |c: &IndelCall| {
        binomial_genotype(c.alt_count as usize, depth as usize, error_rate, tract)
            .is_some_and(|g| !matches!(g.tag, GenotypeTag::HomRef))
    };
    if !stands_alone(first) || !stands_alone(second) {
        return;
    }

    let pooled = (first.alt_count + second.alt_count) as usize;
    if pooled > depth as usize {
        return;
    }
    let (het, hom_alt) = expected_allele_fractions(tract);
    let mass = |p: f64| Binomial::new(depth as usize, clamp_probability(p)).mass(pooled);
    let (p_het, p_compound, p_hom_alt) =
        (mass(het), mass(interpolate(COMPOUND_FRACTION_BY_TRACT, tract)), mass(hom_alt));
    // All three masses can underflow to zero at high depth, and then no `<`
    // comparison holds and every locus would promote to 1/2 with QUAL 0.
    let total = p_het + p_compound + p_hom_alt;
    if total == 0.0 {
        return;
    }
    if p_compound < p_het || p_compound < p_hom_alt {
        return;
    }

    let quality = phred_from_confidence(p_compound / total);
    let (one, two) = (NonZeroU8::new(1).expect("nonzero"), NonZeroU8::new(2).expect("nonzero"));
    calls.truncate(2);
    for call in calls.iter_mut() {
        call.genotype = GenotypeTag::alt_het(one, two);
        call.quality = quality;
    }
}

/// `Binomial::new` panics outside the open unit interval.
fn clamp_probability(p: f64) -> f64 {
    const BOUND: f64 = 1e-9;
    if p.is_finite() { p.clamp(BOUND, 1.0 - BOUND) } else { 0.5 }
}

fn phred_from_confidence(confidence: f64) -> Phred {
    let phred = (-10.0 * (1.0 - confidence).max(1e-300).log10()).min(999.0);
    Phred::from_phred(phred.round().clamp(0.0, 255.0) as u8)
}

struct BinomialGenotype {
    tag: GenotypeTag,
    quality: Phred,
}

/// Classify a site as hom-ref, het, or hom-alt using three binomial hypotheses.
///
/// The het and hom-alt fractions come from [`ALLELE_FRACTION_BY_TRACT`] rather than
/// from 0.5 and `1 - error_rate`: measured, real indels present nowhere near those.
///
/// Returns `None` when depth is zero (no data to genotype).
fn binomial_genotype(
    alt_count: usize,
    total_depth: usize,
    error_rate: f64,
    tract: u32,
) -> Option<BinomialGenotype> {
    if total_depth == 0 {
        return None;
    }

    let alt_one = NonZeroU8::new(1).expect("1 is non-zero");
    let (het, hom_alt) = expected_allele_fractions(tract);

    let p_hom_ref = Binomial::new(total_depth, clamp_probability(error_rate)).mass(alt_count);
    let p_het = Binomial::new(total_depth, clamp_probability(het)).mass(alt_count);
    let p_hom_alt = Binomial::new(total_depth, clamp_probability(hom_alt)).mass(alt_count);

    let total = p_hom_ref + p_het + p_hom_alt;
    if total == 0.0 {
        return None;
    }

    let (best_p, tag) = if p_hom_ref >= p_het && p_hom_ref >= p_hom_alt {
        (p_hom_ref, GenotypeTag::hom_ref())
    } else if p_het >= p_hom_alt {
        (p_het, GenotypeTag::ref_het(alt_one))
    } else {
        (p_hom_alt, GenotypeTag::hom_alt(alt_one))
    };

    let p_best = best_p / total;
    let phred = (-10.0 * (1.0 - p_best).max(1e-300).log10()).min(999.0);

    Some(BinomialGenotype {
        tag,
        quality: Phred::from_phred(phred.round().clamp(0.0, 255.0) as u8),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::pileup::indels::IndelAlleleCounts;
    use seqair_types::Base;

    const ALT_1: NonZeroU8 = NonZeroU8::new(1).expect("nonzero");
    const ALT_2: NonZeroU8 = NonZeroU8::new(2).expect("nonzero");
    const SIMPLE: u32 = 1;
    const LONG_TRACT: u32 = 18;

    fn tag(alt: usize, depth: usize, tract: u32) -> GenotypeTag {
        binomial_genotype(alt, depth, 0.05, tract).expect("non-zero depth").tag
    }

    fn call(seq: &str, alt_count: u32) -> IndelCall {
        IndelCall {
            allele: IndelAllele::Deletion(seq.bytes().map(Base::from).collect()),
            genotype: GenotypeTag::ref_het(ALT_1),
            quality: Phred::from_phred(30_u8),
            ml: None,
            depth: 30,
            alt_count,
            one_sided: false,
        }
    }

    fn counts(ot: u32, ob: u32) -> IndelAlleleCounts {
        IndelAlleleCounts {
            allele: IndelAllele::Deletion([Base::A].into_iter().collect()),
            ot,
            ob,
            unknown_strand: 0,
            noisy: 0,
        }
    }

    /// The fitted curve falls with tract length, and hom-alt falls faster than het —
    /// which is what closes the gap between them inside a repeat.
    #[test]
    fn allele_fractions_fall_with_tract_length() {
        let (het_simple, hom_simple) = expected_allele_fractions(2);
        let (het_long, hom_long) = expected_allele_fractions(18);
        assert!((het_simple - 0.440).abs() < 1e-9 && (hom_simple - 0.900).abs() < 1e-9);
        assert!(het_long < het_simple && hom_long < hom_simple);
        assert!(hom_simple - hom_long > het_simple - het_long);
    }

    /// The failure this exists to fix: real hom-alt indels in a long tract present at
    /// VAF 0.6-0.8, which a hypothesis fixed at `1 - error_rate` calls het.
    #[test]
    fn hom_alt_in_a_long_tract_is_not_called_het() {
        assert_eq!(tag(21, 30, LONG_TRACT), GenotypeTag::HomAlt(ALT_1));
        assert_eq!(tag(19, 30, LONG_TRACT), GenotypeTag::HomAlt(ALT_1));
    }

    /// And the correction must not swallow genuine heterozygotes.
    #[test]
    fn genuine_hets_stay_het() {
        assert_eq!(tag(11, 30, LONG_TRACT), GenotypeTag::RefHet(ALT_1));
        assert_eq!(tag(13, 30, SIMPLE), GenotypeTag::RefHet(ALT_1));
        assert_eq!(tag(0, 30, SIMPLE), GenotypeTag::HomRef);
    }

    #[test]
    fn the_boundary_moves_with_the_tract() {
        let boundary = |tract| {
            (0..=30).find(|&a| tag(a, 30, tract) == GenotypeTag::HomAlt(ALT_1)).unwrap_or(30) as f64
                / 30.0
        };
        assert!(boundary(LONG_TRACT) < boundary(SIMPLE));
    }

    /// Two real alleles at the pooled fraction a 1/2 actually shows (~0.56 in simple
    /// sequence, not ~1.0 — reference bias hits both alleles).
    #[test]
    fn two_real_alleles_become_compound_het() {
        let mut calls = vec![call("A", 9), call("AA", 8)];
        resolve_compound_het(&mut calls, 30, 0.05, SIMPLE);
        assert!(calls.iter().all(|c| c.genotype == GenotypeTag::AltHet(ALT_1, ALT_2)));
    }

    /// A het with a trace second allele pools to about the same fraction, so each
    /// allele has to stand on its own too.
    #[test]
    fn a_het_with_a_trace_second_allele_stays_het() {
        let mut calls = vec![call("A", 14), call("AA", 2)];
        resolve_compound_het(&mut calls, 30, 0.05, SIMPLE);
        assert!(!matches!(calls[0].genotype, GenotypeTag::AltHet(..)));
    }

    /// At high depth all three pooled masses underflow to zero, and then neither
    /// `<` comparison holds. Without a total-mass guard the locus promotes to 1/2
    /// on a hypothesis every model rates as impossible.
    #[test]
    fn an_underflowing_pool_is_not_a_compound_het() {
        let mut calls = vec![call("A", 6927), call("AA", 6927)];
        resolve_compound_het(&mut calls, 20000, 0.05, SIMPLE);
        assert!(calls.iter().all(|c| !matches!(c.genotype, GenotypeTag::AltHet(..))));
    }

    /// The model reinstates a hom-ref allele it rates highly, and its score is
    /// cleared from calls the binomial already accepted so it cannot veto them.
    #[test]
    fn rescue_reinstates_hom_ref_the_model_believes() {
        let probability = |p: f64| Probability::new(p).expect("valid probability");
        let threshold = Some(probability(0.5));
        let mut kept = call("A", 12);
        kept.ml = Some(probability(0.01));
        let mut rescued = call("AA", 3);
        rescued.genotype = GenotypeTag::HomRef;
        rescued.ml = Some(probability(0.9));
        let mut dropped = call("AAA", 2);
        dropped.genotype = GenotypeTag::HomRef;
        dropped.ml = Some(probability(0.1));

        let mut calls = vec![kept, rescued, dropped];
        rescue_hom_ref(&mut calls, threshold);

        assert_eq!(calls.len(), 2, "the unconvincing hom-ref allele is dropped");
        assert_eq!(calls[0].genotype, GenotypeTag::RefHet(ALT_1));
        assert!(calls[0].ml.is_none(), "a low score must not veto a call the binomial kept");
        assert_eq!(calls[1].genotype, GenotypeTag::RefHet(ALT_1), "rescued as a het");
    }

    /// A one-sided allele is emitted for visibility but is not evidence, so it must
    /// not consume a compound-het slot.
    #[test]
    fn a_one_sided_allele_does_not_make_a_compound_het() {
        let mut one_sided = call("AA", 8);
        one_sided.one_sided = true;
        let mut calls = vec![call("A", 9), one_sided];
        resolve_compound_het(&mut calls, 30, 0.05, SIMPLE);
        assert!(calls.iter().all(|c| !matches!(c.genotype, GenotypeTag::AltHet(..))));
    }

    /// A pool that looks hom-alt is one allele with a slippage shadow beside it.
    #[test]
    fn a_hom_alt_pool_is_not_a_compound_het() {
        let mut calls = vec![call("A", 14), call("AA", 13)];
        resolve_compound_het(&mut calls, 30, 0.05, SIMPLE);
        assert!(!matches!(calls[0].genotype, GenotypeTag::AltHet(..)));
    }

    /// Artefacts are strongly enriched among one-sided alleles; measured on TAPS
    /// this removes 14-16 false positives per true call lost.
    #[test]
    fn one_sided_alleles_are_rejected() {
        assert!(!counts(8, 0).on_both_strands());
        assert!(!counts(0, 8).on_both_strands());
        assert!(counts(1, 7).on_both_strands());
    }

    /// Noise comes off both sides of the ratio, so it cannot move the fraction on
    /// its own.
    #[test]
    fn noise_is_excluded_symmetrically() {
        let mut allele = counts(6, 6);
        allele.noisy = 4;
        assert_eq!(allele.clean_total(), 8);
    }
}
