# VCF Fields
Rastair's output follows the [VCFv4.5 specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf).
## Filters
| Name | Description |
| -- | -- |
| **`PASS`** | All filters pass |
| **`lowDp`** | Low read depth |
| **`dnCpG_lowDp`** | Low read depth for de-novo CpG candidate |
| **`dnCpG_bq`** | Low base quality for de-novo CpG candidate |
| **`dnCpG_mapq`** | Low mapping quality for de-novo CpG candidate |
| **`dnCpG_vaf`** | Low variant allele frequency for de-novo CpG candidate |
| **`m_vaf`** | Low variant allele frequency for methylation candidate |
| **`m_bq_ratio`** | Low quality ratio for methylation candidate |
| **`m_pos`** | Alt allele evidence from read edges for methylation candidate |
| **`m_highDp`** | Excessive coverage for methylation candidate |
| **`low_ml_score`** | Machine Learning module prediction below threshold |
## Info Fields
| Name | Description | VCF Type | Rust Type | Occurance |
| -- | -- | -- | -- | -- |
| **`AD`** | Total read depth for each allele | `Integer` | `usize` | R |
| **`BQ`** | RMS base quality | `Float` | `RootMeanSquare` | 1 |
| **`DP`** | Combined depth across samples | `Integer` | `usize` | 1 |
| **`MQ`** | RMS mapping quality | `Float` | `RootMeanSquare` | 1 |
| **`MQ0`** | Number of MAPQ == 0 reads | `Integer` | `usize` | 1 |
| **`NS`** | Number of samples with data | `Integer` | `usize` | 1 |
| **`AS_SB_OT`** | OT counts per allele | `Integer` | `u32` | R |
| **`AS_SB_OB`** | OB counts per allele | `Integer` | `u32` | R |
| **`SC5`** | 5-base sequence context centered on the variant position | `String` | `SmolStr` | 1 |
| **`AF`** | Allele frequency for each ALT allele in the same order as listed (estimated from primary data, not called genotypes) | `Float` | `f64` | A |
| **`ABQ`** | RMS Base quality per allele | `Float` | `f64` | R |
| **`AMQ`** | RMS Map quality per allele | `Float` | `f64` | R |
| **`AS_SS_BQ`** | Strand-specific RMS of base quality per allele (tuples of [reads_ot, reads_ob] for each allele) | `Float` | `f32` | . |
| **`AS_SS_MQ`** | Strand-specific RMS of mapping quality per allele (tuples of [reads_ot, reads_ob] for each allele) | `Float` | `f32` | . |
| **`PIR`** | RMS of relative position in read | `Float` | `f64` | R |
| **`ENT100`** | Shannon entropy of 100bp sequence context around variant position. Value range (0..2) | `Float` | `f64` | 1 |
| **`NAB`** | RMS of number of aligned bases | `Float` | `f64` | R |
| **`NOI`** | RMS of number of indels | `Float` | `f64` | R |
| **`CPG`** | Is this a CpG site? | `Flag` | `bool` | 0 |
| **`CPGnovo`** | De-novo CPG candidate: Could the alt alleles create a new CpG site? | `Flag` | `bool` | 0 |
## Format Fields
| Name | Description | VCF Type | Rust Type | Occurance |
| -- | -- | -- | -- | -- |
| **`GT`** | Genotype | `String` | `GenotypeAllele` | 1 |
| **`GL`** | Genotype likelihoods, Phred-scaled | `Integer` | `Phred>` | G |
| **`GC`** | Genotype confidence, Phred-scaled | `Integer` | `Phred>` | G |
| **`DP`** | Read depth | `Integer` | `usize` | 1 |
| **`M5mC`** | Methylation level at CpG sites | `Float` | `Option<f64>` | 1 |
| **`ML`** | Prediction of methylation/variant likelyhood by Rastair's by machine learning model | `Float` | `f64` | A |
