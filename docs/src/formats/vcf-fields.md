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
| **`dnCpG_adj`** | Included as adjacent position for de-novo CpG candidate, but other position did not pass filters |
| **`m_vaf`** | Low variant allele frequency for methylation candidate |
| **`m_bq_ratio`** | Low quality ratio for methylation candidate |
| **`m_pos`** | Alt allele evidence from read edges for methylation candidate |
| **`m_highDp`** | Excessive coverage for methylation candidate |
| **`pre_ml`** | Low amount of usable evidence, skipping ML |
| **`low_ml_score`** | Machine Learning module prediction below threshold |
## Info Fields
| Name | Description | VCF Type | Rust Type | Occurrence |
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
| **`AS_SS_BQ_OT`** | Strand-specific RMS of base quality per allele on the original top strand | `Float` | `f32` | R |
| **`AS_SS_BQ_OB`** | Strand-specific RMS of base quality per allele on the original bottom strand | `Float` | `f32` | R |
| **`AS_SS_MQ_OT`** | Strand-specific RMS of mapping quality per allele on the original top strand | `Float` | `f32` | R |
| **`AS_SS_MQ_OB`** | Strand-specific RMS of mapping quality per allele on the original bottom strand | `Float` | `f32` | R |
| **`PIR`** | RMS of relative position in read | `Float` | `f64` | R |
| **`ENT100`** | Shannon entropy of 100bp sequence context around variant position. Value range (0..2) | `Float` | `f64` | 1 |
| **`NAB`** | RMS of number of aligned bases | `Float` | `f64` | R |
| **`NOI`** | RMS of number of indels | `Float` | `f64` | R |
| **`M5mC_Strands`** | Number of reads that are evidence for unmodified, modified, no SNP, snp | `Integer` | `u32` | 4 |
| **`CPG`** | Is this a CpG site? | `Flag` | `bool` | 0 |
| **`CPGnovo`** | De-novo CPG candidate: Could the alt alleles create a new CpG site? | `Flag` | `bool` | 0 |
## Format Fields
| Name | Description | VCF Type | Rust Type | Occurrence |
| -- | -- | -- | -- | -- |
| **`GT`** | Genotype | `String` | `GenotypeAllele` | 1 |
| **`GL`** | Genotype likelihoods, Phred-scaled | `Float` | `Phred` | G |
| **`GC`** | Genotype confidence, Phred-scaled | `Float` | `Phred` | G |
| **`DP`** | Read depth | `Integer` | `usize` | 1 |
| **`M5mC`** | Methylation level at CpG sites | `Float` | `Option<f64>` | M |
| **`DPM5mC`** | Total read depth for 5-methylcytosine detection | `Integer` | `Option<u32>` | M |
| **`ADM5mC`** | Read depth supporting 5-methylcytosine modification | `Integer` | `Option<u32>` | M |
| **`ML`** | Prediction of methylation/variant likelihood by rastair's machine learning model | `Float` | `f64` | A |
