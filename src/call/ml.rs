//! Random Forest based classification of positions
//!
//! This module implements feature extraction for machine learning classification
//! of methylation sites in CpG contexts. The feature extraction logic replicates
//! Ben's R notebook analysis for training random forest models on VCF data.
//!
//! ## CpG Feature Extraction
//!
//! The `cpg` module contains functions to extract features from VCF records
//! including:
//! - Basic variant information (ref, alt, mapping quality, etc.)
//! - Sequence context one-hot encoding
//! - Normalized allele depths and strand bias counts
//! - Base and mapping quality metrics by strand
//! - Adjacent position features for methylation evidence
//!
//! The extracted features are returned as an ndarray Array1<f64> suitable
//! for input to a random forest classifier.

#![allow(unused)]

pub use cpg::params_from_record;
pub use models::{MachineLearning, MlResult};

mod models {
    use biosphere::RandomForest;
    use color_eyre::{Result, eyre::Context};
    use std::{fmt, io::Read};
    use tracing::instrument;

    use crate::vcf::Record;

    pub struct MachineLearning {
        threshold: f64,
        cpg: Option<Box<RandomForest>>,
        denovo_cpg: Option<Box<RandomForest>>,
        others: Option<Box<RandomForest>>,
    }

    #[derive(Clone, serde::Serialize, serde::Deserialize)]
    pub enum MlResult {
        None,
        Prediction { prediction: f64, threshold: f64 },
    }

    impl MlResult {
        pub fn pass(&self) -> bool {
            match self {
                MlResult::None => false,
                MlResult::Prediction { prediction, threshold } => prediction >= threshold,
            }
        }
    }

    impl fmt::Debug for MlResult {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            match self {
                MlResult::None => f.debug_tuple("MlResult::None").finish(),
                MlResult::Prediction { prediction, threshold } => f
                    .debug_tuple(if self.pass() { "MlResult::PASS" } else { "MlResult::FAIL" })
                    .field(prediction)
                    .finish(),
            }
        }
    }

    impl std::ops::Deref for MlResult {
        type Target = f64;

        fn deref(&self) -> &Self::Target {
            match self {
                MlResult::None => &0.0,
                MlResult::Prediction { prediction, .. } => prediction,
            }
        }
    }

    impl MachineLearning {
        pub fn disabled() -> Self {
            Self { threshold: 1., cpg: None, denovo_cpg: None, others: None }
        }

        pub fn with_threshold(threshold: f64) -> Self {
            Self {
                threshold,
                cpg: Some(Box::new(
                    load_rf(&include_bytes!("../../models/BS_RF_800-2_CpG.rf.mpk.lz4")[..])
                        .expect("Failed to load CpG RF model"),
                )),
                denovo_cpg: None,
                others: None,
            }
        }

        pub fn cpg(
            &self,
            record: &Record,
            before: Option<&Record>,
            after: Option<&Record>,
        ) -> MlResult {
            let Some(model) = self.cpg.as_ref() else {
                return MlResult::None;
            };
            let features = super::cpg::params_from_record(record, before, after);
            let prediction = model.predict(&features.view().insert_axis(ndarray::Axis(0)));
            match prediction.get(0).copied() {
                Some(p) => MlResult::Prediction { prediction: p, threshold: self.threshold },
                None => MlResult::None,
            }
        }

        pub fn denovo_cpg(
            &self,
            record: &Record,
            before: Option<&Record>,
            after: Option<&Record>,
        ) -> Option<f64> {
            todo!()
        }

        pub fn others(
            &self,
            record: &Record,
            before: Option<&Record>,
            after: Option<&Record>,
        ) -> Option<f64> {
            todo!()
        }
    }

    #[instrument(level = "debug", skip_all)]
    fn load_rf(reader: impl Read) -> Result<RandomForest> {
        let decompress = lz4::Decoder::new(reader).wrap_err("Failed to create LZ4 decoder")?;
        rmp_serde::from_read(decompress).wrap_err("Failed to deserialize random forest")
    }
}

mod cpg {
    use crate::{
        utils::Base,
        vcf::{ByStrand, Record, utils::NoStrandBiasForBaseErrorExt},
    };
    use ndarray::Array1;
    use tracing::{debug, instrument};

    // Ben: the way that these are derived is in the vcf_to_train notebook
    // basically, the `_adj` features are the adjacent G or C (depending on the ref position)
    // p1-p5[ACGT] is just a one-hot encoding of A/C/G/T of the sequence context position 1,2,[ref],4,5
    // all the count features are normalised to depth  or strand-depth
    struct Params {
        ad_alt: f64,
        ad_alt_adj: f64,
        ad_ref: f64,
        alt_a: f64,
        alt_c: f64,
        alt_g: f64,
        alt_t: f64,
        alt_score: f64,
        alt_score_adj: f64,
        bq_alt: f64,
        bq_ob_alt: f64,
        bq_ob_ref: f64,
        bq_ot_alt: f64,
        bq_ot_ref: f64,
        bq_ref: f64,
        is_snp: f64,
        mapq: f64,
        mq_alt: f64,
        mq_ob_alt: f64,
        mq_ob_ref: f64,
        mq_ot_alt: f64,
        mq_ot_ref: f64,
        mq_ref: f64,
        num_aligned_bases_alt: f64,
        num_aligned_bases_ref: f64,
        num_indels_alt: f64,
        num_indels_ref: f64,
        num_mapq0: f64,
        p1_a: f64,
        p1_c: f64,
        p1_g: f64,
        p1_t: f64,
        p2_a: f64,
        p2_c: f64,
        p2_g: f64,
        p2_t: f64,
        p4_a: f64,
        p4_c: f64,
        p4_g: f64,
        p4_t: f64,
        p5_a: f64,
        p5_c: f64,
        p5_g: f64,
        p5_t: f64,
        pos_in_read_alt: f64,
        pos_in_read_ref: f64,
        ref_a: f64,
        ref_c: f64,
        ref_g: f64,
        ref_t: f64,
        region_entropy: f64,
        sb_ob_alt: f64,
        sb_ob_ref: f64,
        sb_ot_alt: f64,
        sb_ot_ref: f64,
    }

    /// Extract feature parameters from a VCF record for CpG classification
    ///
    /// This function implements the feature extraction logic from Ben's R notebook
    /// for training a random forest classifier on CpG methylation sites.
    ///
    /// For CpG methylation analysis, this function specifically looks for:
    /// - C→T transitions (for C reference positions)
    /// - G→A transitions (for G reference positions)
    ///
    /// If multiple alternative alleles are present, only the methylation-relevant
    /// alt (T for C ref, A for G ref) is used for feature extraction.
    ///
    /// # Arguments
    /// * `record` - The current VCF record to extract features from
    /// * `before` - Optional VCF record at position-1 (for adjacent G→A evidence)
    /// * `after` - Optional VCF record at position+1 (for adjacent C→T evidence)
    ///
    /// # Returns
    /// An Array1<f64> containing all the features in the order defined by the Params struct
    #[instrument(level = "debug", skip_all, fields(pos = record.main.pos))]
    pub fn params_from_record(
        record: &Record,
        before: Option<&Record>,
        after: Option<&Record>,
    ) -> Array1<f64> {
        let ref_base = &record.main.r#ref;
        let depth = *record.info.read_depth as f64;

        // Extract basic variant info
        let is_snp = 0.0;

        let mapq = *record.info.mapping_quality;
        let num_mapq0 = *record.info.mapping_quality0 as f64;
        let region_entropy = record.info.entropy[0];

        // Extract sequence context and one-hot encode
        let seq_ctx = &record.info.sequence_context;
        let (p1a, p1c, p1g, p1t) = one_hot_encode_base(seq_ctx.before_2);
        let (p2a, p2c, p2g, p2t) = one_hot_encode_base(seq_ctx.before_1);
        let (p4a, p4c, p4g, p4t) = one_hot_encode_base(seq_ctx.after_1);
        let (p5a, p5c, p5g, p5t) = one_hot_encode_base(seq_ctx.after_2);

        // One-hot encode ref and alt (fetch specific alt for CpG methylation)
        let (ref_a, ref_c, ref_g, ref_t) =
            one_hot_encode_base(Some(Base::from(ref_base.as_bytes()[0])));
        let target_alt = if ref_base == "C" { "T" } else { "A" };
        let (alt_a, alt_c, alt_g, alt_t) =
            if let Some(alt) = record.main.alt.iter().find(|a| a.as_str() == target_alt) {
                one_hot_encode_base(Some(Base::from(alt.as_bytes()[0])))
            } else {
                debug!(target_alt, "No relevant alt allele found for methylation");
                (0.0, 0.0, 0.0, 0.0)
            };

        // Extract normalized allele depths
        let ad_ref = record.info.allele_read_depth.first().copied().unwrap_or(0) as f64 / depth;
        let target_alt_base = if ref_base == "C" { Base::T } else { Base::A };
        let alt_index = record.main.alt.iter().position(|a| a.as_str() == target_alt).unwrap_or(0);
        let ad_alt =
            record.info.allele_read_depth.get(alt_index + 1).copied().unwrap_or(0) as f64 / depth;

        // Extract normalized strand bias counts
        let ref_strand = record.strand_count(Base::from(ref_base.as_bytes()[0])).or_empty();
        let alt_strand = record.strand_count(target_alt_base).or_empty();

        let sb_ot_ref = ref_strand.ot as f64 / depth;
        let sb_ob_ref = ref_strand.ob as f64 / depth;
        let sb_ot_alt = alt_strand.ot as f64 / depth;
        let sb_ob_alt = alt_strand.ob as f64 / depth;

        // Calculate alt_score based on ref base (C vs G)
        let alt_score = if ref_base == "C" {
            // For C: use "ob" (original bottom) strand data
            let bq_ob_alt = get_strand_base_quality(record, target_alt_base).ob;
            let bq_ob_ref = get_strand_base_quality(record, Base::from(ref_base.as_bytes()[0])).ob;
            (sb_ob_alt * bq_ob_alt + 1.0) / (sb_ob_ref * bq_ob_ref + 1.0)
        } else {
            // For G: use "ot" (original top) strand data
            let bq_ot_alt = get_strand_base_quality(record, target_alt_base).ot;
            let bq_ot_ref = get_strand_base_quality(record, Base::from(ref_base.as_bytes()[0])).ot;
            (sb_ot_alt * bq_ot_alt + 1.0) / (sb_ot_ref * bq_ot_ref + 1.0)
        };

        // Extract base quality metrics
        let bq_ref = record.info.allele_base_quality.first().copied().unwrap_or(0.0);
        let bq_alt = record.info.allele_base_quality.get(alt_index + 1).copied().unwrap_or(0.0);
        let bq_ot_ref = get_strand_base_quality(record, Base::from(ref_base.as_bytes()[0])).ot;
        let bq_ob_ref = get_strand_base_quality(record, Base::from(ref_base.as_bytes()[0])).ob;
        let bq_ot_alt = get_strand_base_quality(record, target_alt_base).ot;
        let bq_ob_alt = get_strand_base_quality(record, target_alt_base).ob;

        // Extract mapping quality metrics
        let mq_ref = record.info.allele_map_quality.first().copied().unwrap_or(0.0);
        let mq_alt = record.info.allele_map_quality.get(alt_index + 1).copied().unwrap_or(0.0);
        let mq_ot_ref = get_strand_map_quality(record, Base::from(ref_base.as_bytes()[0])).ot;
        let mq_ob_ref = get_strand_map_quality(record, Base::from(ref_base.as_bytes()[0])).ob;
        let mq_ot_alt = get_strand_map_quality(record, target_alt_base).ot;
        let mq_ob_alt = get_strand_map_quality(record, target_alt_base).ob;

        // Extract other metrics
        let pos_in_read_ref = record.info.position_in_read.first().copied().unwrap_or(0.0);
        let pos_in_read_alt =
            record.info.position_in_read.get(alt_index + 1).copied().unwrap_or(0.0);
        let num_aligned_bases_ref = record.info.num_aligned_bases.first().copied().unwrap_or(0.0);
        let num_aligned_bases_alt =
            record.info.num_aligned_bases.get(alt_index + 1).copied().unwrap_or(0.0);
        let num_indels_ref = record.info.num_indels.first().copied().unwrap_or(0.0);
        let num_indels_alt = record.info.num_indels.get(alt_index + 1).copied().unwrap_or(0.0);

        // Calculate adjacent position features
        let (ad_alt_adj, alt_score_adj) = calculate_adjacent_features(record, before, after);

        let params = Params {
            ad_alt,
            ad_alt_adj,
            ad_ref,
            alt_a,
            alt_c,
            alt_g,
            alt_t,
            alt_score,
            alt_score_adj,
            bq_alt,
            bq_ob_alt,
            bq_ob_ref,
            bq_ot_alt,
            bq_ot_ref,
            bq_ref,
            is_snp,
            mapq,
            mq_alt,
            mq_ob_alt,
            mq_ob_ref,
            mq_ot_alt,
            mq_ot_ref,
            mq_ref,
            num_aligned_bases_alt,
            num_aligned_bases_ref,
            num_indels_alt,
            num_indels_ref,
            num_mapq0,
            p1_a: p1a,
            p1_c: p1c,
            p1_g: p1g,
            p1_t: p1t,
            p2_a: p2a,
            p2_c: p2c,
            p2_g: p2g,
            p2_t: p2t,
            p4_a: p4a,
            p4_c: p4c,
            p4_g: p4g,
            p4_t: p4t,
            p5_a: p5a,
            p5_c: p5c,
            p5_g: p5g,
            p5_t: p5t,
            pos_in_read_alt,
            pos_in_read_ref,
            ref_a,
            ref_c,
            ref_g,
            ref_t,
            region_entropy,
            sb_ob_alt,
            sb_ob_ref,
            sb_ot_alt,
            sb_ot_ref,
        };

        // Convert to Array1<f64>
        params_to_array(params)
    }

    fn one_hot_encode_base(base: Option<Base>) -> (f64, f64, f64, f64) {
        match base {
            Some(Base::A) => (1.0, 0.0, 0.0, 0.0),
            Some(Base::C) => (0.0, 1.0, 0.0, 0.0),
            Some(Base::G) => (0.0, 0.0, 1.0, 0.0),
            Some(Base::T) => (0.0, 0.0, 0.0, 1.0),
            _ => {
                debug!(?base, "Unknown base for one-hot encoding");
                // Unknown or None
                (0.0, 0.0, 0.0, 0.0)
            }
        }
    }

    fn get_strand_base_quality(record: &Record, base: Base) -> ByStrand<f64> {
        record
            .info
            .strand_specific_base_quality
            .iter()
            .find(|x| x.base == base)
            .copied()
            .unwrap_or_default()
    }

    fn get_strand_map_quality(record: &Record, base: Base) -> ByStrand<f64> {
        record
            .info
            .strand_specific_mapping_quality
            .iter()
            .find(|x| x.base == base)
            .copied()
            .unwrap_or_default()
    }

    fn calculate_adjacent_features(
        record: &Record,
        before: Option<&Record>,
        after: Option<&Record>,
    ) -> (f64, f64) {
        let ref_base = &record.main.r#ref;

        if ref_base == "C"
            && let Some(after) = after
            && after.main.r#ref == "G"
            && let Some(alt_index) = after.main.alt.iter().position(|a| a == "A")
        {
            // For C positions: look for G→A transitions in the after record
            let ad_alt =
                after.info.allele_read_depth.get(alt_index + 1).copied().unwrap_or(0) as f64;
            let depth = *after.info.read_depth as f64;
            let ad_alt_norm = ad_alt / depth;

            // Calculate alt_score for G→A
            let alt_strand = after.strand_count(Base::A).or_empty();
            let ref_strand = after.strand_count(Base::G).or_empty();
            let bq_ot_alt = get_strand_base_quality(after, Base::A).ot;
            let bq_ot_ref = get_strand_base_quality(after, Base::G).ot;
            let alt_score =
                (alt_strand.ot as f64 * bq_ot_alt + 1.0) / (ref_strand.ot as f64 * bq_ot_ref + 1.0);

            (ad_alt_norm, alt_score)
        } else if ref_base == "G"
            && let Some(before) = before
            && before.main.r#ref == "C"
            && let Some(alt_index) = before.main.alt.iter().position(|a| a == "T")
        {
            // For G positions: look for C→T transitions in the before record
            let ad_alt =
                before.info.allele_read_depth.get(alt_index + 1).copied().unwrap_or(0) as f64;
            let depth = *before.info.read_depth as f64;
            let ad_alt_norm = ad_alt / depth;

            // Calculate alt_score for C→T
            let alt_strand = before.strand_count(Base::T).or_empty();
            let ref_strand = before.strand_count(Base::C).or_empty();
            let bq_ob_alt = get_strand_base_quality(before, Base::T).ob;
            let bq_ob_ref = get_strand_base_quality(before, Base::C).ob;
            let alt_score =
                (alt_strand.ob as f64 * bq_ob_alt + 1.0) / (ref_strand.ob as f64 * bq_ob_ref + 1.0);

            (ad_alt_norm, alt_score)
        } else {
            // No adjacent evidence for methylation, return defaults
            debug!(%ref_base, before=%before.is_some(), after=%after.is_some(), "No adjacent evidence for methylation");
            (0.0, 0.0)
        }
    }

    fn params_to_array(params: Params) -> Array1<f64> {
        Array1::from(vec![
            params.ad_alt,
            params.ad_alt_adj,
            params.ad_ref,
            params.alt_a,
            params.alt_c,
            params.alt_g,
            params.alt_t,
            params.alt_score,
            params.alt_score_adj,
            params.bq_alt,
            params.bq_ob_alt,
            params.bq_ob_ref,
            params.bq_ot_alt,
            params.bq_ot_ref,
            params.bq_ref,
            params.is_snp,
            params.mapq,
            params.mq_alt,
            params.mq_ob_alt,
            params.mq_ob_ref,
            params.mq_ot_alt,
            params.mq_ot_ref,
            params.mq_ref,
            params.num_aligned_bases_alt,
            params.num_aligned_bases_ref,
            params.num_indels_alt,
            params.num_indels_ref,
            params.num_mapq0,
            params.p1_a,
            params.p1_c,
            params.p1_g,
            params.p1_t,
            params.p2_a,
            params.p2_c,
            params.p2_g,
            params.p2_t,
            params.p4_a,
            params.p4_c,
            params.p4_g,
            params.p4_t,
            params.p5_a,
            params.p5_c,
            params.p5_g,
            params.p5_t,
            params.pos_in_read_alt,
            params.pos_in_read_ref,
            params.ref_a,
            params.ref_c,
            params.ref_g,
            params.ref_t,
            params.region_entropy,
            params.sb_ob_alt,
            params.sb_ob_ref,
            params.sb_ot_alt,
            params.sb_ot_ref,
        ])
    }
}
