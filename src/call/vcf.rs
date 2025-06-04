use rastair2_vcf::{standard_fields::*, *};

// TODO: Ideas for filters
// - from VCF spec
//   filter!(q10, "Quality below 10");
//   filter!(s50, "Less than 50% of samples have data");
// - custom filters
//   filter!(strand_bias, "Significant strand bias detected");
//   filter!(low_coverage, "Low coverage detected");
//   filter!(read_pos, "Variants clustered at read ends");

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
    ],
    // Call data
    //
    // NOTE: The first sub-field must always be the genotype (GT) if it is present.
    format: [SampleReadDepth],
    // hint to allocate this many slots for format data
    min_samples: 1
);
