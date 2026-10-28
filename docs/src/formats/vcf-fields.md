# VCF output fields

## FILTER

| ID | Description |
|----|-------------|
| PASS | All filters passed |
| pre_ml | Low amount of usable evidence, skipping ML |
| m_highDp | Excessive coverage for methylation candidate |
| low_ml_score | Machine Learning module prediction below threshold |
| low_indel_qual | Indel allele the pair-HMM does not call present |
| no_genotype_copy | SNV allele the model passes but the pair-HMM's genotype holds no copy of |

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
| M5mC_Strands | 4 | Integer | Number of reads that are evidence for unmodified, modified (the reads M5mC counts), no SNP, SNP. Always reported; 0,0,0,0 where neither CPG nor CPGnovo is set. |
| CPG | 0 | Flag | Is this a CpG site? |
| CPGnovo | 0 | Flag | De-novo CPG candidate: Could the alt alleles create a new CpG site? |

## FORMAT

| ID | Number | Type | Description |
|----|--------|------|-------------|
| GT | 1 | String | Genotype |
| PS | 1 | Integer | Phase set identifier: the position of the first variant in the phase block, or of its next variant where another block's PS is that position |
| PQ | 1 | Integer | Phasing quality: Phred-scaled probability that the alleles are ordered incorrectly |
| GQ | 1 | Integer | Genotype quality: Phred-scaled chance the genotype is wrong, under a flat prior over the record's genotypes, capped at 99 |
| PL | G | Integer | Phred-scaled genotype likelihoods, the likeliest genotype 0, capped at 255 |
| DP | 1 | Integer | Read depth |
| M5mC | . | Float | Methylation level at CpG sites, one value per CpG context (the reference CpG first), missing where no read shows it |
| DPM5mC | . | Integer | Total read depth for 5-methylcytosine detection, one value per CpG context |
| ADM5mC | . | Integer | Read depth supporting 5-methylcytosine modification, one value per CpG context |
| ML | A | Float | Prediction of methylation/variant likelihood by rastair's machine learning model |
