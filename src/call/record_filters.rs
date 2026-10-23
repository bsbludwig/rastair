use crate::{metrics::PileupMetrics, utils::cli};

/// Filters to apply when deciding whether to output a VCF record
///
/// Here's a table with all the combinations:
///
/// | `vcf_all` | `cpgs_only` | Output behavior                                              |
/// | --------- | ----------- | ------------------------------------------------------------ |
/// | ``        | ``          | All covered CpG sites and variants that PASS                 |
/// | ``        | `-c`        | All covered CpG sites and PASSing de-novo CpGs               |
/// | `--all`   | ``          | All reference CpGs (incl. uncovered) and all variants        |
/// | `--all`   | `-c`        | All reference CpGs (incl. uncovered) and called de-novo CpGs |
///
/// Note: We alwasy report both positions of a CpG.
#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct RecordFilters {
    /// Output all variants, even if they do not pass filters.
    ///
    /// Positions where nothing but the reference was observed are still only
    /// reported when they are a CpG or a de-novo CpG. Every reference CpG is
    /// reported, including those without coverage.
    ///
    /// If combined with `--cpgs-only`, only CpG positions are reported: every
    /// reference CpG including those without coverage, plus de-novo CpGs that
    /// were actually called. De-novo candidates rejected by filters or ML are
    /// not included.
    #[arg(long = "all")]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub vcf_all: bool,

    /// Report CpGs only and default to BED output
    ///
    /// Only report positions that are CpGs in the reference or variants that
    /// would result in a de-novo CpG.
    ///
    /// If combined with `--all`, reference CpGs without coverage in the sample
    /// are also reported. De-novo CpG candidates that were rejected by filters
    /// or ML are not reported.
    #[arg(short = 'c', long, default_value_t = false)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub cpgs_only: bool,
}

/// Everything [`RecordFilters::keeps`] asks about a position.
///
/// A column is answerable before it becomes a [`PileupMetrics`] — which is how
/// the column loop rejects one before finishing it — so the question is
/// asked of these three bits and not of the finished record.
#[derive(Debug, Clone, Copy)]
pub struct PreFilterInputs {
    pub has_alts: bool,
    pub cpg: bool,
    pub has_indels: bool,
}

impl RecordFilters {
    /// Is a position with this evidence worth processing at all?
    pub fn keeps(&self, inputs: PreFilterInputs) -> bool {
        if self.cpgs_only {
            // Filter out pileups that are not CpG if requested
            inputs.cpg
        } else {
            // Otherwise, keep all variant evidence + methylation evidence + indel evidence
            inputs.has_alts || inputs.cpg || inputs.has_indels
        }
    }

    /// Whether every reference CpG is reported, including those no read
    /// reaches, so that two samples' outputs list the same positions.
    pub fn reports_uncovered_cpgs(&self) -> bool {
        self.vcf_all
    }

    /// Filter out pileups that are not relevant based on caller parameters
    ///
    /// Only to speed up processing.
    pub fn pre_filter(&self, pileup: &PileupMetrics) -> bool {
        self.keeps(PreFilterInputs {
            has_alts: !pileup.alts.is_empty(),
            cpg: *pileup.pos_metrics.cpg || pileup.forms_denovo(),
            has_indels: pileup.indel_data.is_some(),
        })
    }

    /// Check if a pileup matches the filter criteria
    pub fn matches(&self, record: &PileupMetrics) -> bool {
        let t = &record.tags;
        debug_assert!(t.set, "record tags not set before filtering");

        let cpg = *record.pos_metrics.cpg || t.denovo_cpg || t.denovo_cpg_partner;

        match (self.vcf_all, self.cpgs_only) {
            (false, false) => t.covered && (cpg || t.variant),
            (false, true) => t.covered && cpg,
            (true, false) => true, // --all means all
            (true, true) => cpg,
        }
    }
}
