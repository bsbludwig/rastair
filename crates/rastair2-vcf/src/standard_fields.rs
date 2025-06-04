//! Standard VCF INFO and FORMAT fields

use crate::{FormatFieldNumber, InfoFieldNumber, format_field, info_field};

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

format_field!(SampleReadDepth(usize), "DP", "Read depth", FormatFieldNumber::Num(1));

mod strand_bias;
pub use strand_bias::StrandBias;
