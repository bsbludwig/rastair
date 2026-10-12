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
| low_indel_qual | Indel allele the pair-HMM does not call present |

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
| PS | 1 | Integer | Phase set identifier: the position of the first variant in the phase block |
| GL | G | Float | Genotype likelihoods, Phred-scaled |
| GC | G | Float | Genotype confidence, Phred-scaled |
| DP | 1 | Integer | Read depth |
| M5mC | . | Float | Methylation level at CpG sites, one value per CpG context |
| DPM5mC | . | Integer | Total read depth for 5-methylcytosine detection, one value per CpG context |
| ADM5mC | . | Integer | Read depth supporting 5-methylcytosine modification, one value per CpG context |
| HPS | 1 | Integer | Phase set the HPC counts are resolved in: either a PS this file emits, whose haplotype 1 it shares, or the position of an unphased heterozygote, in which case haplotype 1 carries the first allele of that site's GT |
| HPC | 4 | Integer | Reads on the TAPS-confounded strand by haplotype: unconverted and converted base on haplotype 1, then on haplotype 2. The confounded pair is C/T at a reference C or T (top strand) and G/A at a reference G or A (bottom strand) |
| HPX | 2 | Integer | Whether each haplotype has a CpG at this position at all: 1 where it carries the cytosine (or the guanine, at the other side of the CpG), 0 where a variant took it away, and 0,0 at a column that is no CpG at all. Omitted where the column could not tell. A haplotype marked 0 converts every read by construction, so its HPC share is its genotype and not a methylation level |
| HPF | 4 | Integer | Fragments of the phase set HPS names, counted once each in HPC's order: unmethylated then methylated on haplotype 1, then the same pair on haplotype 2. Each fragment votes once by the majority of its CpG observations in the set, counted only on a haplotype that carries the CpG at that column, and a tie does not vote. The same four numbers appear on every record of the set |
| ML | A | Float | Prediction of methylation/variant likelihood by rastair's machine learning model |
