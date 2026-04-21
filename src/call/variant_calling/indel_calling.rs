use super::GenotypeTag;
use crate::call::pileup::indels::{IndelAllele, IndelCounts};
use better_default::Default;
use probability::prelude::{Binomial, Discrete as _};
use rastair_types::Phred;
use std::num::NonZeroU8;
use tracing::{instrument, trace};

/// CLI parameters for indel calling.
#[derive(Debug, Clone, clap::Args, serde::Serialize, serde::Deserialize, Default)]
pub struct IndelParams {
    /// Minimum alternate observations to call an indel
    #[arg(long, default_value_t = 2)]
    #[arg(help_heading = crate::utils::cli::sections::FILTER)]
    #[default(2)]
    pub min_indel_ao: u32,

    /// Minimum depth to call an indel
    #[arg(long, default_value_t = 2)]
    #[arg(help_heading = crate::utils::cli::sections::FILTER)]
    #[default(2)]
    pub min_indel_depth: u32,

    /// Error rate for indel genotyping (higher than SNV due to alignment uncertainty)
    #[arg(long, default_value_t = 0.05)]
    #[arg(help_heading = crate::utils::cli::sections::PROCESSING)]
    #[default(0.05)]
    pub indel_error_rate: f64,
}

/// Result of indel calling at a single position for one allele.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IndelCall {
    pub allele: IndelAllele,
    pub genotype: GenotypeTag,
    pub quality: Phred,
    /// Filtered depth (total minus `depth_offset`).
    pub depth: u32,
    /// Reads supporting this indel allele.
    pub alt_count: u32,
}

/// Call indels at a position. Returns empty vec if no indels pass filters.
#[instrument(level = "trace", skip_all)]
pub fn call_indels(indels: &IndelCounts, params: &IndelParams) -> Vec<IndelCall> {
    let mut calls = Vec::new();

    if indels.is_empty() {
        return calls;
    }

    let total_reads = indels.ref_count + indels.total_indel_reads();
    let filtered_depth = total_reads.saturating_sub(indels.depth_offset);

    if filtered_depth < params.min_indel_depth {
        return calls;
    }

    for allele_counts in &indels.alleles {
        // maybe too strict?
        if !allele_counts.on_both_strands() {
            trace!(
                allele = ?allele_counts.allele,
                fwd = allele_counts.fwd,
                rev = allele_counts.rev,
                "Indel skipped: not on both strands"
            );
            continue;
        }

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

        let Some(genotype) =
            binomial_genotype(alt_count as usize, filtered_depth as usize, params.indel_error_rate)
        else {
            continue;
        };

        if matches!(genotype.tag, GenotypeTag::HomRef) {
            trace!(
                allele = ?allele_counts.allele,
                alt_count,
                depth = filtered_depth,
                "Indel skipped: genotyped as hom ref"
            );
            continue;
        }

        calls.push(IndelCall {
            allele: allele_counts.allele.clone(),
            genotype: genotype.tag,
            quality: genotype.quality,
            depth: filtered_depth,
            alt_count,
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

    Some(BinomialGenotype { tag, quality: Phred::from_phred(phred.round() as i32) })
}
