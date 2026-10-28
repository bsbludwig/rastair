//! The options of methylation calling; the calling itself is
//! `metrics::methylation`.

use crate::utils::cli;
use better_default::Default;
use std::num::NonZeroU32;

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct MethylationCallingParams {
    /// Rescue soft-clipped CpG-partner bases
    ///
    /// Aligners often soft-clip a read's fringe base(s) when that mismatches
    /// the reference. But in TAPS, a methylated C reads as T, so a leading C→T
    /// mismatch might get soft-clipped, discarding real methylation evidence
    /// (and same for a trailing G->A). With this flag, a single soft-clipped
    /// base immediately next to an aligned base is "rescued" when it is the
    /// missing partner of a reference CpG (ref C on OT, ref G on OB) and
    /// counted as a normal observation.
    ///
    /// Recovered bases bypass read-end masking (`--nOT`/`--nOB`) since they are
    /// fringe bases by definition, but still go through the base-quality and
    /// mapping-quality filters.
    #[arg(long)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub rescue_soft_clip_cpg: bool,

    /// The most reads that may show a methylation candidate
    ///
    /// A methylation candidate is a C→T alt at a CpG's C or a G→A alt at its
    /// G. One that more reads show is rejected as `m_highDp` whatever the
    /// forests score it: a pileup that deep is a repeat or a collapsed
    /// duplication. Its reads still count as methylation evidence.
    #[arg(long, default_value_t = DEFAULT_M_MAX_COVERAGE)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(DEFAULT_M_MAX_COVERAGE)]
    pub m_max_coverage: NonZeroU32,
}

const DEFAULT_M_MAX_COVERAGE: NonZeroU32 = NonZeroU32::new(1000).expect("1000 is not zero");
