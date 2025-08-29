//! Standard VCF INFO and FORMAT fields

use crate::{FormatFieldNumber, InfoFieldNumber, filter, format_field, info_field};

info_field!(
    AlleleReadDepth(usize),
    "AD",
    "Total read depth for each allele",
    InfoFieldNumber::OnePerAltAndRef
);
info_field!(
    AlleleFrequency(f64),
    "AF",
    "Allele frequency for each ALT allele in the same order as listed (estimated from primary data, not called genotypes)",
    InfoFieldNumber::OnePerAlt
);
info_field!(BaseQuality(RootMeanSquare), "BQ", "RMS base quality", 1);
info_field!(ReadDepth(usize), "DP", "Combined depth across samples", 1);
info_field!(MappingQuality(RootMeanSquare), "MQ", "RMS mapping quality", 1);
info_field!(MappingQuality0(usize), "MQ0", "Number of MAPQ == 0 reads", 1);
info_field!(SamplesWithData(usize), "NS", "Number of samples with data", 1);

format_field!(SampleReadDepth(usize), "DP", "Read depth", 1);
format_field!(
    M5mC(f64),
    "M5mC",
    "Fraction of bases with 5-methylcytosine modification",
    FormatFieldNumber::OnePerPossibleBaseModification
);
format_field!(
    DPM5mC(usize),
    "DPM5mC",
    "Total read depth for 5-methylcytosine detection",
    FormatFieldNumber::OnePerPossibleBaseModification
);
format_field!(
    ADM5mC(usize),
    "ADM5mC",
    "Read depth supporting 5-methylcytosine modification",
    FormatFieldNumber::OnePerPossibleBaseModification
);

filter!(PASS, "All filters pass");

mod strand_bias;
use rastair_types::rms::RootMeanSquare;
pub use strand_bias::StrandBias;
mod genotype;
pub use genotype::{Genotype, GenotypeAllele};
