//! VCF record definition
//!
//! This module defines structures that are used for both representing the data
//! in VCF as well as intermediary data structures in rastair itself. This makes
//! sure we put all data interesting to us also in VCF.
//!
//! See
//! <https://github.com/samtools/hts-specs/blob/0d7f8774658f7cee0a4540b0682174e460726432/VCFv4.5.tex>
//! for the VCF spec.

use rastair2_vcf::{standard_fields::*, *};

mod as_strand_bias;
pub use as_strand_bias::{AlleleSpecificStrandBias, StrandCounts};

mod sequence_context;
pub use sequence_context::SequenceContext;

// TODO: Ideas for filters
// - from VCF spec
//   filter!(q10, "Quality below 10");
//   filter!(s50, "Less than 50% of samples have data");
// - custom filters
//   filter!(strand_bias, "Significant strand bias detected");
//   filter!(low_coverage, "Low coverage detected");
//   filter!(read_pos, "Variants clustered at read ends");

info_field!(
    AllelBaseQuality(f64),
    "ABQ",
    "RMS Base quality per allele",
    InfoFieldNumber::OnePerAltAndRef
);
info_field!(
    AllelMapQuality(f64),
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
info_field!(InCpG, "CPG", "Is this a CpG site?");

format_field!(
    GenotypeLikelihood(Option<f64>),
    "GL",
    "Genotype likelihoods",
    FormatFieldNumber::OnePerGenotype
);
format_field!(
    GenotypeConfidence(Option<f64>),
    "GC",
    "Genotype confidence",
    FormatFieldNumber::OnePerGenotype
);
format_field!(Methylated(Option<f64>), "M5mC", "Methylation level at CpG sites", 1);

vcf_record!(
    // pass or why not
    filters: [PASS],
    // general info
    info: [
        ReadDepthPerAllel,
        BaseQuality,
        ReadDepth,
        MappingQuality,
        MappingQuality0,
        SamplesWithData,
        AlleleSpecificStrandBias,
        SequenceContext,
        AllelFrequency,
        AllelBaseQuality,
        AllelMapQuality,
        PositionInRead,
        Entropy,
        NumAlignedBases,
        NumIndels,
        InCpG,
    ],
    // Call data
    //
    // NOTE: The first sub-field must always be the genotype (GT) if it is present.
    format: [Genotype, GenotypeLikelihood, GenotypeConfidence, SampleReadDepth, Methylated],
    // hint to allocate this many slots for format data
    min_samples: 1
);
