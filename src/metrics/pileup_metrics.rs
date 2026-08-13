use crate::{
    call::{
        pileup::{
            Pileup, SimpleRead,
            indels::{IndelAlleleCounts, IndelCounts},
        },
        variant_calling::{EstimatedGenotype, indel_calling::IndelCall},
    },
    metrics::{MethylationEvidenceStrandInfo, PairedCounts, ReadKey},
    utils::{ByStrand, IntoF64, default, logging::ThisIsABug},
    vcf::{InCpG, Methylated},
};
use better_default::Default;
use color_eyre::{
    Result,
    eyre::{Context, bail},
};
use rastair_vcf::VcfFilter;
use seqair_types::{Base, Probability, RmsAccumulator, RootMeanSquare, SmallVec, SmolStr, Strand};
use std::ops::Deref;
use tracing::trace;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PileupMetrics {
    /// The underlying pileup
    pub pileup: Pileup,
    /// Metrics about the position itself
    pub pos_metrics: PositionMetrics,
    /// Filters that apply to the entire pileup
    pub pos_filters: Filters,
    /// Metrics for the reference allele
    pub ref_metrics: AlleleMetrics,
    /// Metrics and filters for each alternative allele
    pub alts: SmallVec<Alt, 2>,
    /// Counts of (`my_base`, `before_base`) pairs by strand
    pub before_counts: PairedCounts,
    /// Counts of (`my_base`, `after_base`) pairs by strand
    pub after_counts: PairedCounts,
    /// "Tags" for this positions, which will become calls
    pub tags: RecordTags,
    /// Aggregated indel counts at this position.
    #[serde(default)]
    pub indels: IndelCounts,
    /// Called indel variants at this position (populated during `process_region`).
    #[serde(default)]
    pub indel_calls: Vec<IndelCall>,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct RecordTags {
    // The tags have been calculated. Mainly to debug :)
    pub set: bool,
    /// This position has coverage
    pub covered: bool,
    /// This is in a CpG site
    pub cpg: bool,
    /// This is a de-novo CpG site (not the partner)
    pub denovo_cpg: bool,
    /// This is the partner position of a de-novo CpG site
    pub denovo_cpg_partner: bool,
    /// This position is a variant (but not a denovo CpG)
    pub variant: bool,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Alt {
    pub base: Base,
    pub metrics: AlleleMetrics,
    pub filters: AltFilters,
    pub call: AltCall,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum AltCall {
    #[default]
    Uncalled,
    RealVariant,
    MethylationEvidenceOnly {
        for_base: Base,
    },
    ReadError,
}

impl PileupMetrics {
    /// Create new metrics from pileup
    ///
    /// NOTE: The extended metrics in `PositionMetrics` are not set here and
    /// need to be set later using `set_extended_metrics`.
    pub fn new(pileup: Pileup) -> Result<Self> {
        let ref_base = pileup.reference_base;
        let total_reads = pileup.reads.len();

        // Single pass: accumulate per-base metrics, position-level RMS, and discover alleles
        let mut accumulators = PerBaseAccumulators::default();
        let mut pos_baseq = RmsAccumulator::new();
        let mut pos_mapq = RmsAccumulator::new();
        let mut mapq0: u32 = 0;
        let mut alt_bases: SmallVec<Base, 4> = SmallVec::new();
        for read in pileup.reads.iter() {
            let qual_sq = f64::from(read.qual).powi(2);
            let mapq_sq = f64::from(read.mapq).powi(2);
            accumulators.accumulate(read, qual_sq, mapq_sq);
            pos_baseq.add_squared(qual_sq);
            pos_mapq.add_squared(mapq_sq);
            if read.mapq == 0 {
                mapq0 += 1;
            }
            if read.base.known_index().is_some()
                && read.base != ref_base
                && !alt_bases.contains(&read.base)
            {
                alt_bases.push(read.base);
            }
        }

        trace!(pos = pileup.pos, ?ref_base, ?alt_bases, "New pileup");

        let pos_metrics = PositionMetrics::from_pileup(
            &pileup,
            PositionMetricsExt::default(),
            pos_baseq.finish(),
            pos_mapq.finish(),
            mapq0,
        );

        // Reference base can be Unknown (N in FASTA) — no accumulator slot exists for it,
        // and no reads will ever match it, so just use default metrics.
        let ref_metrics = if let Some(acc) = accumulators.take(ref_base) {
            acc.finish(ref_base, total_reads, &pileup)
                .wrap_err("Failed to compute allele metrics for reference")?
        } else {
            AlleleMetrics { base: ref_base, ..default() }
        };

        let alts = alt_bases
            .iter()
            .map(|&base| {
                // alt_bases only contains known bases (filtered above), so take always succeeds
                let acc = accumulators
                    .take(base)
                    .ok_or_else(|| color_eyre::eyre::eyre!("unknown base {base} in alt_bases"))?;
                let metrics = acc
                    .finish(base, total_reads, &pileup)
                    .wrap_err("Failed to compute allele metrics for alt")?;
                Ok(Alt { base, metrics, filters: AltFilters::default(), call: default() })
            })
            .collect::<Result<_>>()?;

        let indels = aggregate_indels(&pileup);

        let mut before_counts = PairedCounts::default();
        let mut after_counts = PairedCounts::default();
        for read in pileup.reads.iter() {
            if read.strand == Strand::Unknown {
                continue;
            }
            if let Some(before) = read.before_base {
                before_counts.increment(ReadKey {
                    strand: read.strand,
                    current: read.base,
                    adj: before,
                });
            }
            if let Some(after) = read.after_base {
                after_counts.increment(ReadKey {
                    strand: read.strand,
                    current: read.base,
                    adj: after,
                });
            }
        }

        Ok(PileupMetrics {
            pileup,
            pos_metrics,
            pos_filters: Filters::default(),
            ref_metrics,
            alts,
            before_counts,
            after_counts,
            tags: RecordTags::default(),
            indels,
            indel_calls: Vec::new(),
        })
    }

    /// Get reference base
    pub fn ref_base(&self) -> Base {
        self.pileup.reference_base
    }

    pub fn contig(&self) -> SmolStr {
        self.pileup.region.contig.clone()
    }

    pub fn pos(&self) -> u32 {
        self.pileup.pos
    }

    pub fn contig_pos(&self) -> SmolStr {
        use std::fmt::Write as _;
        let mut res = seqair_types::smol_str::SmolStrBuilder::new();

        write!(&mut res, "{}:{}", self.contig(), self.pos()).expect("works");

        res.finish()
    }

    pub fn alt(&self, alt: Base) -> Option<&AlleleMetrics> {
        self.alts.iter().find(|a| a.base == alt).map(|a| &a.metrics)
    }

    pub fn allele(&self, base: Base) -> Option<&AlleleMetrics> {
        if base == self.ref_base() {
            Some(&self.ref_metrics)
        } else {
            self.alts.iter().find(|a| a.base == base).map(|a| &a.metrics)
        }
    }

    pub fn alt_metrics(&self, alt: Base) -> Option<MetricsForAlt<'_>> {
        let alt = self.alts.iter().find(|a| a.base == alt);
        alt.map(|alt| MetricsForAlt { metrics: self, alt: &alt.metrics })
    }

    pub fn alt_filters(&self, alt: Base) -> Option<&AltFilters> {
        self.alts.iter().find(|a| a.base == alt).map(|a| &a.filters)
    }

    pub fn alt_filters_mut(&mut self, alt: Base) -> Option<&mut AltFilters> {
        self.alts.iter_mut().find(|a| a.base == alt).map(|a| &mut a.filters)
    }

    /// Get all alternative bases in the pileup (for lookup with mutation)
    pub fn alts(&self) -> SmallVec<Base, 4> {
        self.alts.iter().map(|a| a.base).collect()
    }

    pub fn alts_metrics(&self) -> impl Iterator<Item = &AlleleMetrics> {
        self.alts.iter().map(|a| &a.metrics)
    }

    pub fn ref_alts_metrics(&self) -> impl Iterator<Item = &AlleleMetrics> {
        std::iter::once(&self.ref_metrics).chain(self.alts.iter().map(|a| &a.metrics))
    }

    pub fn forms_denovo(&self) -> bool {
        *self.pos_metrics.denovo_adj || self.alts.iter().any(|a| *a.metrics.denovo)
    }

    pub fn pass(&self, ml_threshold: Option<Probability>) -> bool {
        if self.pos_filters.other_pos_in_denovo_passes {
            return true;
        }
        self.pos_filters.pass() && self.alts.iter().any(|a| a.filters.pass(ml_threshold))
    }
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[cfg_attr(test, derive(Default))]
pub struct PositionMetrics {
    /// Read depth, i.e., number of reads covering this position
    pub depth: u32,
    /// Base quality
    pub baseq: RootMeanSquare,
    /// Mapping quality
    pub mapq: RootMeanSquare,
    /// Number of reads with mapping quality 0
    pub mapq0: u32,
    /// Is this position in a CpG context in the reference?
    pub cpg: InCpG,

    /// Extended metrics
    // set by `call` later since they depend on more context
    #[serde(flatten)]
    pub extended: PositionMetricsExt,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct PositionMetricsExt {
    /// Entropy of the surrounding region
    pub region_entropy: f64,
    /// Estimated genotype
    pub genotype: Option<EstimatedGenotype>,
    /// Methylation strand info
    pub methylation_strand_info: MethylationEvidenceStrandInfo,
    /// Methylation beta
    pub methylated: Methylated,
    /// Is this position a de-novo cpg candidate?
    pub denovo_adj: DenovoAdjecent,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum DenovoAdjecent {
    #[default]
    No,
    ThisIsTheMatchingC,
    ThisIsTheMatchingG,
}

impl Deref for DenovoAdjecent {
    type Target = bool;

    fn deref(&self) -> &Self::Target {
        match self {
            DenovoAdjecent::No => &false,
            _ => &true,
        }
    }
}

impl PositionMetrics {
    pub fn from_pileup(
        pileup: &Pileup,
        extended: PositionMetricsExt,
        baseq: RootMeanSquare,
        mapq: RootMeanSquare,
        mapq0: u32,
    ) -> Self {
        PositionMetrics {
            depth: u32::try_from(pileup.reads.len()).expect("depth fits into u32"),
            baseq,
            mapq,
            mapq0,
            cpg: InCpG::from(pileup),

            // These fields are given by all
            extended,
        }
    }
}

impl Deref for PositionMetrics {
    type Target = PositionMetricsExt;

    fn deref(&self) -> &Self::Target {
        &self.extended
    }
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct AlleleMetrics {
    pub base: Base,
    /// Read depth, i.e. number of reads supporting this allele
    pub depth: u32,
    /// base quality for the allele
    pub baseq: RootMeanSquare,
    /// mapping quality for the allele
    pub mapq: RootMeanSquare,
    /// count of reads by strand, also known as strand bias
    pub strand_count: ByStrand<u32>,
    /// base quality by strand
    pub baseq_s: ByStrand<RootMeanSquare>,
    /// mapping quality by strand
    pub mapq_s: ByStrand<RootMeanSquare>,
    /// number of aligned bases in read
    pub num_aligned_bases: RootMeanSquare,
    /// number of indels in read
    pub num_indels: RootMeanSquare,
    /// relative position in read
    pub position_in_read: RootMeanSquare,
    /// Allele frequency
    pub allele_frequency: Probability,
    /// does this alt form a de-novo cpg?
    pub denovo: FormsDenovo,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
pub enum FormsDenovo {
    #[default]
    No,
    ThisBecomesC,
    ThisBecomesG,
}

impl Deref for FormsDenovo {
    type Target = bool;

    fn deref(&self) -> &Self::Target {
        match self {
            FormsDenovo::No => &false,
            _ => &true,
        }
    }
}

impl FormsDenovo {
    pub fn some(&self) -> Option<Self> {
        match self {
            FormsDenovo::No => None,
            _ => Some(*self),
        }
    }
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct AltFilters {
    /// ML prediction: probability this is a true variant
    pub ml: Option<Probability>,
    pub filters: Filters,
}

impl AltFilters {
    pub fn pass(&self, ml_threshold: Option<Probability>) -> bool {
        if self.filters.other_pos_in_denovo_passes {
            return true;
        }
        if let Some(ml_threshold) = ml_threshold
            && let Some(ml) = self.ml
        {
            ml >= ml_threshold
        } else {
            self.filters.is_empty()
        }
    }
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct Filters {
    pub other_pos_in_denovo_passes: bool,
    filters: SmallVec<SmolStr, 6>,
}

impl Filters {
    pub fn add(&mut self, filter: impl VcfFilter, condition: impl FnOnce() -> bool) {
        if condition() {
            self.filters.push(filter.filter());
        }
    }

    pub fn merge(&mut self, other: Filters) {
        for filter in other.filters {
            if !self.filters.contains(&filter) {
                self.filters.push(filter);
            }
        }
    }

    pub fn pass(&self) -> bool {
        self.other_pos_in_denovo_passes || self.filters.is_empty()
    }
}

impl Deref for Filters {
    type Target = SmallVec<SmolStr, 6>;

    fn deref(&self) -> &Self::Target {
        &self.filters
    }
}

#[derive(Debug, Clone, Default)]
struct AlleleAccumulator {
    depth: u32,
    baseq: RmsAccumulator,
    mapq: RmsAccumulator,
    baseq_ot: RmsAccumulator,
    baseq_ob: RmsAccumulator,
    mapq_ot: RmsAccumulator,
    mapq_ob: RmsAccumulator,
    aligned: RmsAccumulator,
    indels: RmsAccumulator,
    pos_in_read: RmsAccumulator,
    ot_count: u32,
    ob_count: u32,
}

impl AlleleAccumulator {
    fn add(&mut self, read: &SimpleRead, qual_sq: f64, mapq_sq: f64) {
        self.depth += 1;
        self.baseq.add_squared(qual_sq);
        self.mapq.add_squared(mapq_sq);
        match read.strand {
            Strand::OT => {
                self.ot_count += 1;
                self.baseq_ot.add_squared(qual_sq);
                self.mapq_ot.add_squared(mapq_sq);
            }
            Strand::OB => {
                self.ob_count += 1;
                self.baseq_ob.add_squared(qual_sq);
                self.mapq_ob.add_squared(mapq_sq);
            }
            Strand::Unknown => {}
        }
        self.aligned.add(f64::from(read.matching_bases));
        self.indels.add(f64::from(read.indels));
        self.pos_in_read.add(f64::from(read.position.pos) / f64::from(read.position.read_length));
    }

    fn finish(self, base: Base, total_reads: usize, pileup: &Pileup) -> Result<AlleleMetrics> {
        use Base::*;

        if self.depth == 0 {
            // Can happen at canonical CpG sites with no evidence for a particular allele
            trace!(
                pos = pileup.pos,
                ref_base = ?pileup.reference_base,
                ?base,
                pileup_reads = pileup.reads.len(),
                "No reads for allele"
            );
            return Ok(AlleleMetrics { base, ..default() });
        }

        // Should be impossible: depth > 0 implies total_reads > 0 since reads were counted
        // from the same pileup. Guard defensively since a division by zero here would produce
        // NaN/Inf which Probability::new rejects anyway, but better to fail with a clear message.
        if total_reads == 0 {
            bail!("allele has depth {} but pileup has 0 total reads — this is a bug", self.depth);
        }

        let denovo = if base == pileup.reference_base {
            FormsDenovo::No
        } else if pileup.ref_before() == C && base == G {
            FormsDenovo::ThisBecomesG
        } else if pileup.ref_after() == G && base == C {
            FormsDenovo::ThisBecomesC
        } else {
            FormsDenovo::No
        };

        Ok(AlleleMetrics {
            base,
            depth: self.depth,
            baseq: self.baseq.finish(),
            mapq: self.mapq.finish(),
            strand_count: ByStrand { base, ot: self.ot_count, ob: self.ob_count },
            baseq_s: ByStrand { base, ot: self.baseq_ot.finish(), ob: self.baseq_ob.finish() },
            mapq_s: ByStrand { base, ot: self.mapq_ot.finish(), ob: self.mapq_ob.finish() },
            num_aligned_bases: self.aligned.finish(),
            num_indels: self.indels.finish(),
            position_in_read: self.pos_in_read.finish(),
            allele_frequency: Probability::new(self.depth.f() / total_reads.f())
                .wrap_err("allele frequency not in [0,1]")
                .this_is_a_bug()?,
            denovo,
        })
    }
}

/// Per-base accumulators indexed by [`Base::known_index`], one slot per `Base::KNOWN`.
#[derive(Debug, Default)]
struct PerBaseAccumulators([AlleleAccumulator; 4]);

impl PerBaseAccumulators {
    fn accumulate(&mut self, read: &SimpleRead, qual_sq: f64, mapq_sq: f64) {
        let Some(idx) = read.base.known_index() else { return };
        self.0[idx].add(read, qual_sq, mapq_sq);
    }

    fn take(&mut self, base: Base) -> Option<AlleleAccumulator> {
        let idx = base.known_index()?;
        Some(std::mem::take(&mut self.0[idx]))
    }
}

pub struct MetricsForAlt<'p> {
    pub metrics: &'p PileupMetrics,
    pub alt: &'p AlleleMetrics,
}

impl MetricsForAlt<'_> {
    pub fn is_evidence_for_methylation(&self) -> bool {
        (self.metrics.pos_metrics.cpg == InCpG::C && self.alt.base == Base::T)
            || (self.metrics.pos_metrics.cpg == InCpG::G && self.alt.base == Base::A)
    }
}

pub struct MetricsForIndel<'p> {
    pub metrics: &'p PileupMetrics,
    pub indel: &'p IndelCall,
}

fn aggregate_indels(pileup: &Pileup) -> IndelCounts {
    if pileup.indel_observations.is_empty() {
        return IndelCounts {
            ref_count: pileup.reads.len() as u32,
            noisy_ref_count: pileup.noisy_ref_count,
            ..Default::default()
        };
    }

    let mut alleles: SmallVec<IndelAlleleCounts, 2> = SmallVec::new();

    for obs in &pileup.indel_observations {
        match alleles.iter_mut().find(|e| e.allele == obs.allele) {
            Some(entry) => entry.add(obs.strand, obs.noisy),
            None => {
                let mut entry = IndelAlleleCounts::new(obs.allele.clone());
                entry.add(obs.strand, obs.noisy);
                alleles.push(entry);
            }
        }
    }

    let total_indel_reads: u32 = alleles.iter().map(|a| a.total()).sum();
    let ref_count = (pileup.reads.len() as u32).saturating_sub(total_indel_reads);

    // Strand mix of the locus as a whole, used as the null for the per-allele
    // strand-bias test. Counted over `reads` — the same post-dedup fragments the
    // observations came from — so the two are on one granularity.
    let mut ot_depth = 0;
    let mut ob_depth = 0;
    for read in pileup.reads.iter() {
        match read.strand {
            Strand::OT => ot_depth += 1,
            Strand::OB => ob_depth += 1,
            Strand::Unknown => {}
        }
    }

    IndelCounts {
        alleles,
        ot_depth,
        ob_depth,
        ref_count,
        depth_offset: pileup.depth_offset,
        noisy_ref_count: pileup.noisy_ref_count,
    }
}

#[cfg(test)]
mod aggregate_indel_tests {
    use super::*;
    use crate::call::pileup::SimpleReads;
    use crate::call::pileup::indels::{IndelAllele, IndelObservation};
    use crate::sequence::{ChunkRegion, Region};
    use crate::utils::{SequenceContext, default};

    fn observation(strand: Strand, reverse: bool) -> IndelObservation {
        IndelObservation {
            allele: IndelAllele::Insertion([Base::A].into_iter().collect()),
            strand,
            reverse,
            pos_in_read: 10,
            read_length: 100,
            mapq: 60,
            base_qual: 40,
            matching_bases: 100,
            num_indels_in_read: 1,
            insertion_base_quals: default(),
            post_del_base_qual: 0,
            has_repeat: false,
            noisy: false,
        }
    }

    fn pileup_with(observations: Vec<IndelObservation>) -> Pileup {
        Pileup {
            region: ChunkRegion {
                region: Region { contig: "chrT".into(), start: 1000, end: 1100 },
                last_position: 2000,
                overlap_start: 0,
                overlap_end: 0,
            },
            context: SequenceContext::default(),
            pos: 1002,
            reads: SimpleReads(Vec::new().into()),
            reference_base: Base::T,
            indel_observations: observations.into_iter().collect(),
            depth_offset: 0,
            homopolymer_run: 0,
            dinucleotide_run: 0,
            soft_clip_count: 0,
            noisy_ref_count: 0,
            indel_ref_window: default(),
            indel_ref_anchor: 0,
        }
    }

    /// Regression guard: per-fragment deduplication keeps one mate per fragment, and
    /// for a coordinate-sorted FR library that is always the forward one — so every
    /// surviving observation has `reverse == false`. The strand split must still
    /// show both strands, which it only does if it splits on OT/OB rather than on
    /// the alignment reverse flag.
    #[test]
    fn strand_split_survives_fragment_dedup() {
        let counts = aggregate_indels(&pileup_with(vec![
            observation(Strand::OT, false),
            observation(Strand::OT, false),
            observation(Strand::OB, false),
        ]));

        let allele = counts.alleles.first().expect("one allele");
        assert_eq!((allele.ot, allele.ob), (2, 1));
        assert_eq!(allele.total(), 3);
        assert_eq!(
            allele.strand_bias_p_value(0.5),
            1.0,
            "a 2/1 split is as balanced as an odd count allows"
        );
    }

    /// The mirror of the above: the reverse flag carries no strand information
    /// after deduplication, so two fragments from the same duplex strand must not
    /// look balanced merely because their reverse flags differ.
    #[test]
    fn opposite_reverse_flags_on_one_strand_are_not_a_split() {
        let counts = aggregate_indels(&pileup_with(vec![
            observation(Strand::OT, false),
            observation(Strand::OT, true),
        ]));

        let allele = counts.alleles.first().expect("one allele");
        assert_eq!((allele.ot, allele.ob), (2, 0));
    }

    /// Reads whose orientation could not be determined still count toward depth
    /// but contribute nothing to the strand split either way.
    #[test]
    fn unknown_strand_counts_toward_total_only() {
        let counts = aggregate_indels(&pileup_with(vec![
            observation(Strand::OT, false),
            observation(Strand::Unknown, false),
        ]));

        let allele = counts.alleles.first().expect("one allele");
        assert_eq!((allele.ot, allele.ob, allele.unknown_strand), (1, 0, 1));
        assert_eq!(allele.total(), 2);
    }

    /// The null for the strand-bias test comes from the fragments *not* supporting
    /// the allele, so an allele cannot dilute its own null.
    #[test]
    fn null_strand_fraction_excludes_the_alleles_own_support() {
        let mut counts = aggregate_indels(&pileup_with(vec![
            observation(Strand::OT, false),
            observation(Strand::OT, false),
        ]));
        counts.ot_depth = 12;
        counts.ob_depth = 8;

        let allele = counts.alleles.first().expect("one allele");
        // 10 of the 18 non-supporting fragments are OT, not 12 of 20 — plus one
        // fragment of prior mass centred on the locus' own 12.5/21 OT share.
        assert_eq!(counts.null_ot_fraction(allele), (10.0 + 12.5 / 21.0) / 19.0);
    }

    /// Noise on the supporting fragments has to reach the allele it supports, or
    /// the hard-filter pathway would drop noisy fragments from the denominator
    /// while leaving them in the numerator — the asymmetry the split was meant to
    /// remove.
    #[test]
    fn supporting_fragment_noise_reaches_its_allele() {
        let mut noisy = observation(Strand::OB, false);
        noisy.noisy = true;
        let mut pileup = pileup_with(vec![observation(Strand::OT, false), noisy]);
        pileup.noisy_ref_count = 3;
        // 2 supporting + 8 non-supporting fragments, 3 of the latter noisy.
        pileup.reads = SimpleReads(vec![SimpleRead::default(); 10].into());

        let counts = aggregate_indels(&pileup);
        let allele = counts.alleles.first().expect("one allele");
        assert_eq!((allele.total(), allele.noisy, allele.clean_total()), (2, 1, 1));
        assert_eq!(counts.ref_count, 8);
        assert_eq!(counts.noisy_ref_count, 3);
        assert_eq!(counts.clean_depth(), 6, "5 clean reference + 1 clean supporting");
    }

    /// A homozygous indel leaves nothing to estimate the null from; the test must
    /// fall back to unbiased rather than divide by zero.
    #[test]
    fn null_strand_fraction_falls_back_when_every_fragment_supports_the_allele() {
        let mut counts = aggregate_indels(&pileup_with(vec![
            observation(Strand::OT, false),
            observation(Strand::OT, false),
        ]));
        // A locus covered on OT alone, every fragment of it supporting the allele:
        // the homozygous case. With no background left, the null is the locus' own
        // mix, so the allele is not judged against a strand it was never read on.
        counts.ot_depth = 2;
        counts.ob_depth = 0;

        let allele = counts.alleles.first().expect("one allele");
        assert!(
            counts.null_ot_fraction(allele) > 0.8,
            "null must follow the locus, not collapse to a balanced 0.5"
        );
        assert!(allele.strand_bias_p_value(counts.null_ot_fraction(allele)) > 0.05);
    }
}
