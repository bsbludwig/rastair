# VCF output fields

## FILTER

| ID | Description |
|----|-------------|
| PASS | All filters passed |
| lowDp | Low read depth |
| dnCpG_lowDp | Low read depth for de-novo CpG candidate |
| dnCpG_bq | Low base quality for de-novo CpG candidate |
| dnCpG_mapq | Low mapping quality for de-novo CpG candidate |
| dnCpG_vaf | Low variant allele frequency for de-novo CpG candidate |
| dnCpG_adj | Included as adjacent position for de-novo CpG candidate, but other position did not pass filters |
| m_vaf | Low variant allele frequency for methylation candidate |
| m_bq_ratio | Low quality ratio for methylation candidate |
| m_pos | Alt allele evidence from read edges for methylation candidate |
| m_highDp | Excessive coverage for methylation candidate |
| pre_ml | Low amount of usable evidence, skipping ML |
| low_ml_score | Machine Learning module prediction below threshold |
| indel_strand | Indel allele supported on only one bisulfite strand |

## INFO

| ID | Number | Type | Description |
|----|--------|------|-------------|
| AD | R | Integer | Total read depth for each allele |
| BQ | 1 | Float | RMS base quality |
| DP | 1 | Integer | Combined depth across samples |
| MQ | 1 | Float | RMS mapping quality |
| MQ0 | 1 | Integer | Number of MAPQ == 0 reads |
| NS | 1 | Integer | Number of samples with data |
| AS_SB_OT | R | Integer | OT counts per allele |
| AS_SB_OB | R | Integer | OB counts per allele |
| SC5 | 1 | String | 5-base sequence context centered on the variant position |
| AF | A | Float | Allele frequency for each ALT allele in the same order as listed (estimated from primary data, not called genotypes) |
| ABQ | R | Float | RMS Base quality per allele |
| AMQ | R | Float | RMS Map quality per allele |
| AS_SS_BQ_OT | R | Float | Strand-specific RMS of base quality per allele on the original top strand |
| AS_SS_BQ_OB | R | Float | Strand-specific RMS of base quality per allele on the original bottom strand |
| AS_SS_MQ_OT | R | Float | Strand-specific RMS of mapping quality per allele on the original top strand |
| AS_SS_MQ_OB | R | Float | Strand-specific RMS of mapping quality per allele on the original bottom strand |
| PIR | R | Float | RMS of relative position in read |
| ENT100 | 1 | Float | Shannon entropy of 100bp sequence context around variant position. Value range (0..2) |
| NAB | R | Float | RMS of number of aligned bases |
| NOI | R | Float | RMS of number of indels |
| M5mC_Strands | 4 | Integer | Number of reads that are evidence for unmodified, modified, no SNP, SNP. Always reported, can be non-zero while CPG and CPGnovo are unset. |
| CPG | 0 | Flag | Is this a CpG site? |
| CPGnovo | 0 | Flag | De-novo CPG candidate: Could the alt alleles create a new CpG site? |

## FORMAT

| ID | Number | Type | Description |
|----|--------|------|-------------|
| GT | 1 | String | Genotype |
| GL | G | Float | Genotype likelihoods, Phred-scaled |
| GC | G | Float | Genotype confidence, Phred-scaled |
| DP | 1 | Integer | Read depth |
| M5mC | . | Float | Methylation level at CpG sites, one value per CpG context |
| DPM5mC | . | Integer | Total read depth for 5-methylcytosine detection, one value per CpG context |
| ADM5mC | . | Integer | Read depth supporting 5-methylcytosine modification, one value per CpG context |
| ML | A | Float | Prediction of methylation/variant likelihood by rastair's machine learning model |
