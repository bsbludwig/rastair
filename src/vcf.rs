//! VCF record definition
//!
//! This module defines structures that are used for both representing the data
//! in VCF as well as intermediary data structures in rastair itself. This makes
//! sure we put all data interesting to us also in VCF.
//!
//! See
//! <https://github.com/samtools/hts-specs/blob/0d7f8774658f7cee0a4540b0682174e460726432/VCFv4.5.tex>
//! for the VCF spec.

use rastair_vcf::{standard_fields::*, *};
use seqair_types::Phred;

mod as_strand_bias;
pub use as_strand_bias::AlleleSpecificStrandBias;
mod asq;
pub use asq::{StrandSpecificBaseQuality, StrandSpecificMappingQuality};
mod cpg;
pub use cpg::InCpG;
mod denovo_cpg;
pub use denovo_cpg::DeNovoCpGCandidate;
mod methylation;
pub use methylation::{CpgBeta, CpgOrigin, Methylated, MethylationAltDepth, MethylationDepth};

mod metrics_to_vcf;

use crate::metrics::MethylationEvidenceStrandInfo;
pub use crate::utils::SequenceContext;

filter!(lowDp, "Low read depth");
filter!(dnCpG_lowDp, "Low read depth for de-novo CpG candidate");
filter!(dnCpG_bq, "Low base quality for de-novo CpG candidate");
filter!(dnCpG_mapq, "Low mapping quality for de-novo CpG candidate");
filter!(dnCpG_vaf, "Low variant allele frequency for de-novo CpG candidate");
filter!(
    dnCpG_adj,
    "Included as adjacent position for de-novo CpG candidate, but other position did not pass filters"
);
filter!(m_vaf, "Low variant allele frequency for methylation candidate");
filter!(m_bq_ratio, "Low quality ratio for methylation candidate");
filter!(m_pos, "Alt allele evidence from read edges for methylation candidate");
filter!(m_highDp, "Excessive coverage for methylation candidate");
filter!(low_ml_score, "Machine Learning module prediction below threshold");
filter!(pre_ml, "Low amount of usable evidence, skipping ML");
filter!(indel_strand, "Indel allele not supported on both strands");
filter!(indel_hom_ref, "Indel genotyped as homozygous reference by binomial model");

info_field!(
    AlleleBaseQuality(f64),
    "ABQ",
    "RMS Base quality per allele",
    InfoFieldNumber::OnePerAltAndRef
);
info_field!(
    AlleleMapQuality(f64),
    "AMQ",
    "RMS Map quality per allele",
    InfoFieldNumber::OnePerAltAndRef
);
info_field!(
    PositionInRead(f64),
    "PIR",
    "RMS of relative position in read",
    InfoFieldNumber::OnePerAltAndRef
);
info_field!(
    Entropy(f64),
    "ENT100",
    "Shannon entropy of 100bp sequence context around variant position. Value range (0..2)",
    InfoFieldNumber::Num(1)
);
info_field!(
    NumAlignedBases(f64),
    "NAB",
    "RMS of number of aligned bases",
    InfoFieldNumber::OnePerAltAndRef
);
info_field!(NumIndels(f64), "NOI", "RMS of number of indels", InfoFieldNumber::OnePerAltAndRef);

format_field!(
    GenotypeLikelihood(Option<Phred>),
    "GL",
    "Genotype likelihoods, Phred-scaled",
    FormatFieldNumber::OnePerGenotype
);
format_field!(
    GenotypeConfidence(Option<Phred>),
    "GC",
    "Genotype confidence, Phred-scaled",
    FormatFieldNumber::OnePerGenotype
);
format_field!(
    MachineLearningPrediction(f64),
    "ML",
    "Prediction of methylation/variant likelihood by rastair's machine learning model",
    FormatFieldNumber::OnePerAlt
);

vcf_record!(
    // pass or why not
    filters: [
        PASS,
        lowDp,
        dnCpG_lowDp, dnCpG_bq, dnCpG_mapq, dnCpG_vaf, dnCpG_adj,
        m_vaf, m_bq_ratio, m_pos, m_highDp,
        pre_ml, low_ml_score,
        indel_strand, indel_hom_ref,
    ],
    // general info
    // Fields marked "default" are included in minimal VCF output by default
    // Others can be enabled via --vcf-info-fields CLI flag
    info: [
        AlleleReadDepth default,
        BaseQuality default,
        ReadDepth default,
        MappingQuality default,
        MappingQuality0,
        SamplesWithData,
        AlleleSpecificStrandBias,
        SequenceContext,
        AlleleFrequency,
        AlleleBaseQuality,
        AlleleMapQuality,
        StrandSpecificBaseQuality,
        StrandSpecificMappingQuality,
        PositionInRead,
        Entropy,
        NumAlignedBases,
        NumIndels,
        MethylationEvidenceStrandInfo default,
        InCpG default,
        DeNovoCpGCandidate default,
    ],
    // Call data
    //
    // NOTE: The first sub-field must always be the genotype (GT) if it is present.
    // Fields marked "default" are included in minimal VCF output by default
    // Others can be enabled via --vcf-format-fields CLI flag
    format: [
        Genotype default,
        GenotypeLikelihood default,
        GenotypeConfidence default,
        SampleReadDepth default,
        Methylated default,
        MethylationDepth default,
        MethylationAltDepth default,
        MachineLearningPrediction default
    ],
    // hint to allocate this many slots for format data
    min_samples: 1
);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_field_config_with_field_ids() {
        let info_ids: Vec<InfoFieldId> = ["AF", "MQ0"].iter().map(|s| s.parse().unwrap()).collect();
        let config = FieldConfig::default().with_field_ids(&info_ids, &[]);

        // Default fields should still be enabled
        assert!(config.info.read_depth, "DP is default");
        assert!(config.info.allele_read_depth, "AD is default");
        assert!(config.format.genotype, "GT is default");
        assert!(config.format.genotype_likelihood, "GL is default");
        assert!(config.format.genotype_confidence, "GC is default");
        assert!(config.format.machine_learning_prediction, "ML is default");

        // Additional INFO fields should be enabled
        assert!(config.info.allele_frequency, "AF should be enabled");
        assert!(config.info.mapping_quality0, "MQ0 should be enabled");

        // Other fields should remain disabled
        assert!(!config.info.samples_with_data, "NS should not be enabled");
        assert!(!config.info.allele_base_quality, "ABQ should not be enabled");
    }

    #[test]
    fn test_field_config_invalid_field_id() {
        let result = "INVALID".parse::<InfoFieldId>();
        assert!(result.is_err(), "Should error on invalid INFO field");

        let result = "INVALID".parse::<FormatFieldId>();
        assert!(result.is_err(), "Should error on invalid FORMAT field");
    }

    #[test]
    fn test_field_filtering_integration() {
        use rastair_vcf::{Compression, Contig, VcfBuilder, VcfFormat};
        use tempfile::TempDir;

        // Create a test VCF writer with default config
        let temp_dir = TempDir::new().expect("create temp dir");
        let temp_file = temp_dir.path().join("test_default.vcf");

        let writer = VcfBuilder::new(&temp_file, VcfFormat::Vcf, Compression::Off, 1)
            .expect("create builder");

        let contigs = [Contig { name: seqair_types::SmolStr::new("chr1"), length: 1000 }];
        let samples = [seqair_types::SmolStr::new("sample1")];

        let mut vcf = writer.build::<Record>(&contigs, &samples).expect("build vcf");

        // Check default fields are enabled
        let config = vcf.config_mut();
        assert!(config.info.read_depth, "DP should be in default");
        assert!(config.info.allele_read_depth, "AD should be in default");
        assert!(config.info.base_quality, "BQ should be in default");
        assert!(config.info.mapping_quality, "MQ should be in default");
        assert!(config.info.methylation_evidence_strand_info, "MESI should be in default");
        assert!(config.info.in_cp_g, "InCpG should be in default");
        assert!(config.info.de_novo_cp_g_candidate, "DeNovoCpGCandidate should be in default");

        // Check non-default fields are disabled
        assert!(!config.info.allele_frequency, "AF should not be in default");
        assert!(!config.info.samples_with_data, "NS should not be in default");

        // Create a custom config with additional fields
        let temp_file2 = temp_dir.path().join("test_custom.vcf");
        let writer2 = VcfBuilder::new(&temp_file2, VcfFormat::Vcf, Compression::Off, 1)
            .expect("create builder");

        let info_ids: Vec<InfoFieldId> = ["AF", "MQ0"].iter().map(|s| s.parse().unwrap()).collect();
        let custom_config = FieldConfig::default().with_field_ids(&info_ids, &[]);

        let mut vcf2 = writer2
            .build::<Record>(&contigs, &samples)
            .expect("build vcf")
            .with_config(custom_config);

        // Custom config should have additional fields enabled
        assert!(vcf2.config_mut().info.allele_frequency, "AF should be enabled");
        assert!(vcf2.config_mut().info.mapping_quality0, "MQ0 should be enabled");
        // Default fields should still be enabled
        assert!(vcf2.config_mut().info.read_depth, "DP should still be enabled");
    }
}
