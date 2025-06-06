use rastair2_vcf::{standard_fields::*, *};
use smol_str::SmolStr;

// TODO: Ideas for filters
// - from VCF spec
//   filter!(q10, "Quality below 10");
//   filter!(s50, "Less than 50% of samples have data");
// - custom filters
//   filter!(strand_bias, "Significant strand bias detected");
//   filter!(low_coverage, "Low coverage detected");
//   filter!(read_pos, "Variants clustered at read ends");

info_field!(
    SequenceContext(SmolStr),
    "SC5",
    "5-base sequence context centered on the variant position",
    InfoFieldNumber::Num(1)
);
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
    "Shannon entropy of 100bp sequence context around variant position",
    InfoFieldNumber::Num(1)
);
info_field!(
    NumAlignedBases(f64),
    "NAB",
    "RMS of number of aligned bases",
    InfoFieldNumber::OnePerAltAndRef
);
info_field!(NumIndels(f64), "NOI", "RMS of number of indels", InfoFieldNumber::OnePerAltAndRef);

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
        StrandBias,
        SequenceContext,
        AllelFrequency,
        AllelBaseQuality,
        AllelMapQuality,
        PositionInRead,
        Entropy,
        NumAlignedBases,
        NumIndels,
    ],
    // Call data
    //
    // NOTE: The first sub-field must always be the genotype (GT) if it is present.
    format: [SampleReadDepth],
    // hint to allocate this many slots for format data
    min_samples: 1
);
