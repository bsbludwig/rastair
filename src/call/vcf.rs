use rastair2_vcf::*;

info_field!(BaseQuality(f64), "BQ", "RMS base quality", FieldNumber::Num(1));
info_field!(MappingQuality(f64), "MQ", "RMS mapping quality", FieldNumber::Num(1));

format_field!(ReadDepth(usize), "DP", "Read depth", FieldNumber::Num(1));

// filter!(q10, "Quality below 10");
// filter!(s50, "Less than 50% of samples have data");

vcf_record!(
    filters: [],
    info: [BaseQuality, MappingQuality],
    format: [ReadDepth],
    min_samples: 1
);
