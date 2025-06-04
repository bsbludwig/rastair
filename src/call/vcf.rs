use rastair2_vcf::*;

info_field!(BaseQuality(f64), "BQ", "RMS base quality", FieldNumber::Num(1));
info_field!(MappingQuality(f64), "MQ", "RMS mapping quality", FieldNumber::Num(1));

format_field!(ReadDepth(usize), "DP", "Read depth", FieldNumber::Num(1));

// filter!(q10, "Quality below 10");
// filter!(s50, "Less than 50% of samples have data");

vcf_record!(
    // pass or why not
    filters: [],
    // general info
    info: [BaseQuality, MappingQuality],
    // Call data
    //
    // NOTE: The first sub-field must always be the genotype (GT) if it is present.
    format: [ReadDepth],
    // hint to allocate this many slots for format data
    min_samples: 1
);
