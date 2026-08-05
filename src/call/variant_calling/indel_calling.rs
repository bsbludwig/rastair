use super::GenotypeTag;
use crate::call::pileup::indels::{IndelAllele, IndelCounts};
use better_default::Default;
use probability::prelude::{Binomial, Discrete as _};
use seqair_types::{Phred, Probability};
use std::num::NonZeroU8;
use tracing::{instrument, trace};

// Self-contained non-ML hard-filter indel pathway (used under `--no-ml`).
pub mod hard_filters;

/// CLI parameters for indel calling.
#[derive(Debug, Clone, clap::Args, serde::Serialize, serde::Deserialize, Default)]
pub struct IndelParams {
    /// Enable experimental indel calling
    ///
    /// When disabled, Rastair calls SNPs and methylation only.
    #[arg(long, default_value_t = false)]
    pub experimental_indels: bool,

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
    /// Filtered depth (total minus `depth_offset`).
    pub depth: u32,
    /// Reads supporting this indel allele.
    pub alt_count: u32,
    /// Non-ML hard-filter verdict; `None` on the ML path.
    pub hard_filter_verdict: Option<hard_filters::IndelVerdict>,
}

/// Call indels at a position. Returns empty vec if no indels pass filters.
///
/// When `ml_enabled` is true, the binomial genotype test is used only for
/// informational genotyping — all alleles passing min AO and min depth
/// are forwarded to the ML model. When ML is off, the binomial test
/// acts as a hard gate (`hom_ref` alleles are rejected).
#[instrument(level = "trace", skip_all)]
pub fn call_indels(indels: &IndelCounts, params: &IndelParams, ml_enabled: bool) -> Vec<IndelCall> {
    let mut calls = Vec::new();

    if indels.is_empty() {
        return calls;
    }

    // With ML off, indels go through the non-ML hard-filter chain instead.
    if !ml_enabled {
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

        let genotype =
            binomial_genotype(alt_count as usize, filtered_depth as usize, params.indel_error_rate);

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
/// Returns `None` when depth is zero (no data to genotype).
fn binomial_genotype(
    alt_count: usize,
    total_depth: usize,
    error_rate: f64,
) -> Option<BinomialGenotype> {
    if total_depth == 0 {
        return None;
    }

    let alt_one = NonZeroU8::new(1).expect("1 is non-zero");

    let p_hom_ref = Binomial::new(total_depth, error_rate).mass(alt_count);
    let p_het = Binomial::new(total_depth, 0.5).mass(alt_count);
    let p_hom_alt = Binomial::new(total_depth, 1.0 - error_rate).mass(alt_count);

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
