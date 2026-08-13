use super::GenotypeTag;
use crate::call::pileup::indels::{IndelAllele, IndelCounts};
use better_default::Default;
use probability::prelude::{Binomial, Discrete as _};
use seqair_types::{Phred, Probability};
use std::num::NonZeroU8;
use tracing::{instrument, trace};

// Self-contained non-ML hard-filter indel pathway (used by `--experimental-indels`).
pub mod hard_filters;

/// How noisy fragments — soft-clipped, or carrying a terminal tandem repeat —
/// are kept out of the hard-filter pathway's genotype ratio.
///
/// The three settings differ along two independent axes: whether a noisy
/// *supporting* fragment still counts as an observation for the `min_indel_ao`
/// gate, and whether it counts in the numerator of the alt/depth ratio. The
/// denominator drops noisy non-supporting fragments in every mode.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum, serde::Serialize, serde::Deserialize,
)]
#[serde(rename_all = "kebab-case")]
pub enum IndelNoiseExclusion {
    /// Drop noisy fragments from the alternate count and the depth alike.
    ///
    /// Noise is a property of the read, so it lands on supporting and
    /// non-supporting fragments at the same rate; removing it from both sides
    /// leaves VAF where it was. It also removes them from the `min_indel_ao`
    /// gate, so an allele whose support is mostly noisy is never a candidate.
    Symmetric,
    /// Keep noisy fragments as observations for the `min_indel_ao` gate, but
    /// exclude them from both sides of the ratio.
    ///
    /// [`Self::Symmetric`] charges noise twice — once at candidate generation and
    /// again in the ratio — which costs candidates outright rather than merely
    /// genotyping them conservatively. This charges it once.
    RatioOnly,
    /// Drop noisy fragments from the depth only.
    ///
    /// Inflates VAF by roughly the local noise rate, which walks low-VAF alleles
    /// across the binomial's hom-ref/het boundary at VAF ≈ 0.218. Statistically
    /// this is a bias, not a correction — but the inflation is largest in repeat
    /// context, which is also where alignment under-reports genuine indel support,
    /// so it acts as a crude reference-bias correction and measures materially
    /// more sensitive on GIAB. Kept available rather than deleted for that reason.
    DepthOnly,
}

impl IndelNoiseExclusion {
    /// Support counted toward the `min_indel_ao` gate, and toward the numerator of
    /// the genotype ratio, for one allele.
    fn alt_counts(self, allele: &crate::call::pileup::indels::IndelAlleleCounts) -> (u32, u32) {
        match self {
            Self::Symmetric => (allele.clean_total(), allele.clean_total()),
            Self::RatioOnly => (allele.total(), allele.clean_total()),
            Self::DepthOnly => (allele.total(), allele.total()),
        }
    }

    /// Denominator of the genotype ratio at a locus.
    fn depth(self, indels: &IndelCounts) -> u32 {
        match self {
            Self::Symmetric | Self::RatioOnly => indels.clean_depth(),
            Self::DepthOnly => (indels.ref_count + indels.total_indel_reads())
                .saturating_sub(indels.noisy_ref_count),
        }
    }
}

/// CLI parameters for indel calling.
#[derive(Debug, Clone, clap::Args, serde::Serialize, serde::Deserialize, Default)]
pub struct IndelParams {
    /// Enable experimental indel calling using hard filters
    ///
    /// Candidate indels are accepted by a fixed filter chain (depth, alternate
    /// observations, strand bias, binomial genotype) with no ML scoring.
    /// When disabled, Rastair calls SNPs and methylation only.
    #[arg(long, default_value_t = false)]
    pub experimental_indels: bool,

    /// Enable experimental indel calling scored by the ML model
    ///
    /// Selects the machine-learning indel pathway instead of the hard-filter
    /// chain of `--experimental-indels`. Has no effect together with `--no-ml`,
    /// which falls back to the hard-filter chain.
    #[arg(long, default_value_t = false, conflicts_with = "experimental_indels")]
    pub experimental_indels_ml: bool,

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

    /// Expected alternate fraction of a heterozygous indel.
    ///
    /// The 0.5 a heterozygote would show without reference bias. Reads carrying an
    /// indel are harder to place than reads matching the reference, so a genuine
    /// het indel is systematically observed below 0.5 and the binomial calls it
    /// homozygous reference. Lowering this models that bias where it happens.
    ///
    /// It is not interchangeable with `--indel-error-rate`, which moves the same
    /// hom-ref/het boundary from the other side: lowering the error rate weakens
    /// the hom-ref hypothesis against genuine sequencing noise as well, so it buys
    /// recall by giving up precision at noise loci. Lowering this leaves the
    /// hom-ref null where the data says it is.
    #[arg(long, default_value_t = 0.5)]
    #[arg(help_heading = crate::utils::cli::sections::PROCESSING)]
    #[default(0.5)]
    pub indel_het_vaf: f64,

    /// Significance level for the indel strand-bias filter. Off by default.
    ///
    /// An indel allele is rejected when the OT/OB split of its supporting
    /// fragments is this unlikely or less under the strand mix of the rest of the
    /// locus. Lower is more permissive; `0`, the default, disables the filter.
    ///
    /// It defaults to off because the hypothesis it tests is false for TAPS: OT
    /// and OB reads present different sequence after C→T conversion, so genuine
    /// indel support is strand-asymmetric for reasons that are not artifacts. At
    /// the 0.05 it used to default to, measured on chr12 against GIAB HG001, it
    /// rejected 4,288 true indels to remove 252 false ones and cost 0.098 F1. The
    /// p-value is still informative — it belongs in the ML feature vector rather
    /// than in a hard gate.
    #[arg(long, default_value_t = 0.0)]
    #[arg(help_heading = crate::utils::cli::sections::FILTER)]
    #[default(0.0)]
    pub indel_strand_bias_alpha: f64,

    /// How noisy fragments are excluded from the hard-filter genotype ratio
    ///
    /// "Noisy" means soft-clipped or carrying a terminal tandem repeat.
    /// `symmetric` drops them from the alternate count and the depth alike, which
    /// leaves VAF unbiased but also costs candidates at the min-AO gate;
    /// `ratio-only` still counts them as observations for that gate; `depth-only`
    /// drops them from the depth alone, which raises VAF by roughly the local
    /// noise rate and so acts as a sensitivity boost. See
    /// [`IndelNoiseExclusion`](crate::call::variant_calling::indel_calling::IndelNoiseExclusion).
    #[arg(long, value_enum, default_value_t = IndelNoiseExclusion::Symmetric)]
    #[arg(help_heading = crate::utils::cli::sections::FILTER)]
    #[default(IndelNoiseExclusion::Symmetric)]
    pub indel_noise_exclusion: IndelNoiseExclusion,

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

impl IndelParams {
    /// Whether indels are called at all, by either pathway.
    pub fn enabled(&self) -> bool {
        self.experimental_indels || self.experimental_indels_ml
    }

    /// Whether the ML pathway is selected *and* the ML model is available.
    ///
    /// `--experimental-indels-ml --no-ml` degrades to the hard-filter chain
    /// rather than erroring, so the two flags compose the way the rest of the
    /// caller's `--no-ml` handling does.
    pub fn use_ml(&self, ml_enabled: bool) -> bool {
        self.experimental_indels_ml && ml_enabled
    }
}

/// Result of indel calling at a single position for one allele.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IndelCall {
    pub allele: IndelAllele,
    pub genotype: GenotypeTag,
    pub quality: Phred,
    /// ML score: probability this indel is a true variant.
    pub ml: Option<Probability>,
    /// Genotyping depth. The ML path subtracts `depth_offset` from total reads;
    /// the hard-filter path uses `IndelCounts::clean_depth`, which drops noisy
    /// fragments from both sides of the ratio.
    pub depth: u32,
    /// Fragments supporting this indel allele, on the same footing as
    /// [`Self::depth`] — noise-excluded on the hard-filter path, raw on the ML one.
    pub alt_count: u32,
    /// Non-ML hard-filter verdict; `None` on the ML path.
    pub hard_filter_verdict: Option<hard_filters::IndelVerdict>,
}

/// Call indels at a position. Returns empty vec if no indels pass filters.
///
/// On the ML path (`use_ml`) the binomial genotype test is only informational —
/// every allele passing min AO and min depth is forwarded to the ML model. On the
/// hard-filter path each genotyped allele is emitted with an
/// [`hard_filters::IndelVerdict`] recording whether it passed; failing alleles are
/// tagged with a VCF FILTER rather than dropped.
#[instrument(level = "trace", skip_all)]
pub fn call_indels(indels: &IndelCounts, params: &IndelParams, use_ml: bool) -> Vec<IndelCall> {
    let mut calls = Vec::new();

    if indels.is_empty() {
        return calls;
    }

    if !use_ml {
        return hard_filters::call_indels(indels, params);
    }

    let total_reads = indels.ref_count + indels.total_indel_reads();
    let filtered_depth = total_reads.saturating_sub(indels.depth_offset);

    if filtered_depth < params.min_indel_depth {
        return calls;
    }

    for allele_counts in &indels.alleles {
        let alt_count = allele_counts.total();
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
            params.indel_het_vaf,
        );

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
            hard_filter_verdict: None,
        });
    }

    calls
}

struct BinomialGenotype {
    tag: GenotypeTag,
    quality: Phred,
}

/// Classify a site as hom-ref, het, or hom-alt using three binomial hypotheses.
///
/// The heterozygous hypothesis sits at `het_vaf` rather than a fixed 0.5 — see
/// [`IndelParams::indel_het_vaf`]. The homozygous-alternate one is left at
/// `1 - error_rate`: reference bias moves the het/hom-ref boundary, which is where
/// recall is lost, and not the het/hom-alt one.
///
/// Returns `None` when depth is zero (no data to genotype).
fn binomial_genotype(
    alt_count: usize,
    total_depth: usize,
    error_rate: f64,
    het_vaf: f64,
) -> Option<BinomialGenotype> {
    if total_depth == 0 {
        return None;
    }

    let alt_one = NonZeroU8::new(1).expect("1 is non-zero");

    // `Binomial::new` panics outside the open unit interval, and these arrive from
    // the command line.
    const BOUND: f64 = 1e-9;
    let clamp = |p: f64| if p.is_finite() { p.clamp(BOUND, 1.0 - BOUND) } else { 0.5 };

    let p_hom_ref = Binomial::new(total_depth, clamp(error_rate)).mass(alt_count);
    let p_het = Binomial::new(total_depth, clamp(het_vaf)).mass(alt_count);
    let p_hom_alt = Binomial::new(total_depth, clamp(1.0 - error_rate)).mass(alt_count);

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
