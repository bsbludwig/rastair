//! Standard VCF INFO and FORMAT fields

use crate::{FormatFieldNumber, InfoFieldNumber, filter, format_field, info_field};

info_field!(
    ReadDepthPerAllel(usize),
    "AD",
    "Total read depth for each allele",
    InfoFieldNumber::OnePerAltAndRef
);
info_field!(
    AllelFrequency(f64),
    "AF",
    "Allele frequency for each ALT allele in the same order as listed (estimated from primary data, not called genotypes)",
    InfoFieldNumber::OnePerAlt
);
info_field!(BaseQuality(f64), "BQ", "RMS base quality", InfoFieldNumber::Num(1));
info_field!(ReadDepth(usize), "DP", "Combined depth across samples", InfoFieldNumber::Num(1));
info_field!(MappingQuality(f64), "MQ", "RMS mapping quality", InfoFieldNumber::Num(1));
info_field!(MappingQuality0(usize), "MQ0", "Number of MAPQ == 0 reads", InfoFieldNumber::Num(1));
info_field!(SamplesWithData(usize), "NS", "Number of samples with data", InfoFieldNumber::Num(1));

format_field!(SampleReadDepth(usize), "DP", "Read depth", FormatFieldNumber::Num(1));

filter!(PASS, "All filters pass");

mod strand_bias;
pub use strand_bias::StrandBias;
