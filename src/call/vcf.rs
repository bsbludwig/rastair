use rastair2_vcf::{standard_fields::StrandBias, *};

info_field!(
    ReadDepthPerAllel(usize),
    "AD",
    "Total read depth for each allele",
    InfoFieldNumber::OnePerAltAndRef
);
info_field!(BaseQuality(f64), "BQ", "RMS base quality", InfoFieldNumber::Num(1));
info_field!(ReadDepth(usize), "DP", "Combined depth across samples", InfoFieldNumber::Num(1));
info_field!(MappingQuality(f64), "MQ", "RMS mapping quality", InfoFieldNumber::Num(1));
info_field!(MappingQuality0(usize), "MQ0", "Number of MAPQ == 0 reads", InfoFieldNumber::Num(1));
info_field!(SamplesWithData(usize), "NS", "Number of samples with data", InfoFieldNumber::Num(1));
// info_field!(StrandBias(usize), "SB", "Strand bias", InfoFieldNumber::Num(4));

format_field!(SampleReadDepth(usize), "DP", "Read depth", FormatFieldNumber::Num(1));

// filter!(q10, "Quality below 10");
// filter!(s50, "Less than 50% of samples have data");

vcf_record!(
    // pass or why not
    filters: [],
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
