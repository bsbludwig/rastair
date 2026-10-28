use crate::{call::process::UncoveredCpgs, metrics::PileupMetrics, utils::cli};

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
/// Both positions of a CpG are always reported.
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

/// What makes a position worth processing; a column with none of these is
/// dropped.
///
/// The same for every output: CpG-only output drops its records only in
/// [`RecordFilters::matches`], because a non-CpG column still decides CpGs.
/// Two SNVs can form one (`TA>CG`), which only their calls reveal, and the
/// pair-HMM's verdict at a CpG depends on the variants around it.
///
/// A column is answerable before it becomes a [`PileupMetrics`], which is how
/// the column loop rejects one before finishing it.
#[derive(Debug, enumset::EnumSetType)]
pub enum ColumnEvidence {
    Alt,
    Cpg,
    Indel,
}

/// Drop the columns without [`ColumnEvidence`], as the column loop
/// does before finishing them, for the hand-built pileups of the tests:
/// `columns` are every column of a region, so a column's genomic neighbours
/// are still beside it.
#[cfg(test)]
pub fn pre_filter(columns: &mut Vec<PileupMetrics>) {
    use crate::{
        utils::adjacent,
        vcf::{CpgSide, TapsContext},
    };

    let alt_sides = |column: Option<&PileupMetrics>| -> enumset::EnumSet<CpgSide> {
        let inserted = column.iter().flat_map(|c| c.indel_data.iter());
        let inserted = inserted.flat_map(|d| d.counts.alleles.iter().map(|a| &a.allele));
        let alts = column.iter().flat_map(|c| c.alts.iter().map(|alt| alt.base));
        alts.chain(inserted.filter_map(|a| a.last_inserted()))
            .filter_map(CpgSide::of_base)
            .collect()
    };
    let keep: Vec<bool> = columns
        .iter()
        .enumerate()
        .map(|(i, column)| {
            let before = i.checked_sub(1).and_then(|i| columns.get(i));
            let after = columns.get(i + 1);
            let before = alt_sides(before.filter(|b| adjacent(b, column)));
            let after = alt_sides(after.filter(|a| adjacent(column, a)));
            let mut evidence = enumset::EnumSet::<ColumnEvidence>::empty();
            if !column.alts.is_empty() {
                evidence |= ColumnEvidence::Alt;
            }
            if CpgSide::taken_by(
                column.reference_base,
                TapsContext::of(column).sides() | CpgSide::beside(before, after),
            )
            .is_some()
            {
                evidence |= ColumnEvidence::Cpg;
            }
            if column.indel_data.is_some() {
                evidence |= ColumnEvidence::Indel;
            }
            !evidence.is_empty()
        })
        .collect();
    let mut keep = keep.into_iter();
    columns.retain(|_| keep.next().unwrap_or(false));
}

impl RecordFilters {
    /// Whether every reference CpG is reported, including those no read
    /// reaches.
    pub fn uncovered_cpgs(&self) -> UncoveredCpgs {
        if self.vcf_all { UncoveredCpgs::Fill } else { UncoveredCpgs::Skip }
    }

    /// Check if a pileup matches the filter criteria
    pub fn matches(&self, record: &PileupMetrics) -> bool {
        let cpg = record.in_cpg();

        match (self.vcf_all, self.cpgs_only) {
            (false, false) => record.covered() && (cpg || record.is_variant()),
            (false, true) => record.covered() && cpg,
            (true, false) => true, // --all means all
            (true, true) => cpg,
        }
    }
}
