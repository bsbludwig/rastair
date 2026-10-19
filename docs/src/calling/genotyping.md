# Genotyping

Rastair uses @ML scores to determine which @altAllele:pl represent true genetic @variant:pl, then applies statistical methods to estimate zygosity.
Rastair only supports diploid genotyping (two chromosome copies) at this time.

## Multi-Allelic Sites

When multiple alternate alleles are present at the same genomic position, they are combined into a single VCF record with comma-separated ALT values.
For example, with `REF=A` and `ALT=T,G`, the genotype represents the complete set of observed alleles:

- A genotype of `1/2` means "on one chromosome there was a T (allele 1), on the other was a G (allele 2)"
- A genotype of `0/1` means "on one chromosome there was the reference A (allele 0), on the other was a T (allele 1)"
- A genotype of `1/1` means "both chromosomes had a T (allele 1)"

This representation follows VCF specification standards where allele indices start at 0 for the reference, and 1, 2, 3... for the alternate alleles in order.

## Genotype Calls

`0/0` denotes two reference copies, `0/1` one reference and one alternate copy,
`1/1` two copies of the same alternate, and `1/2` two different alternates.
Passing a variant filter makes an allele eligible for genotyping; it does not determine its copy number.

Models trained with `--ml-features haplotypes` automatically enable molecular
genotyping for SNPs. No `--phase` flag is needed for this step. Models using
the Standard feature set retain count-based genotyping. Indels are genotyped
by the [pair-HMM indel caller](indels.md) from its own haplotypes.

For molecular genotyping, PairHMM scores each captured fragment against candidate
chromosome sequences. A diploid pair mixes the two chromosome likelihoods with
equal weight, and independent fragments contribute to the joint likelihood.
Overlapping mates represent one molecule: disagreement reduces the consensus
quality rather than contributing two independent votes.

The genotype posterior combines these likelihoods with the SNP prior and
marginalizes over retained neighboring haplotypes. SNPs consider all ten
unordered A/C/G/T genotypes over the observed alleles. Connected candidates
share a bounded regional search (32 sites, 2 kb, 64 retained haplotypes, 32
buffered alternatives). One molecular likelihood matrix supplies all site
marginals in the region. The best pair for every local genotype and the
all-reference chromosome are protected during pruning. If those protections
exceed the budget, inference splits into separate regions without removing the
earlier candidate's alternatives.

The forest remains an allele filter. Its score is not multiplied into the
molecular genotype posterior, since the forest already sees molecular features.
Rejected alleles cannot appear in the emitted genotype, but their posterior
mass remains uncertainty when calculating confidence. Genotyping runs before
methylation calling and before the final genotype-conditioned phaser.

When molecular evidence is unavailable, Rastair falls back to count-based
genotyping. Examples include masked intervals, indels over 100 bp, unsafe
complex repeat events, wholly aligned-fallback evidence, and competing-indel
sets that exceed the search limit. Molecular scoring uses the existing read
quality and candidate filters; it does not assemble arbitrary complex alleles.

## Confidence Scoring

On SNP records with a genotype estimate and on indel records, FORMAT `GTM`
identifies the method: `pairhmm` or `counts`. FORMAT `GQ` reports model-derived
genotype confidence, capped at 99. If no SNP genotype estimate is available,
these two annotations are omitted; the existing output fallback remains `0/0`.

The forest decides whether an allele is present and the molecular likelihood
decides its dosage. A molecular genotype is chosen among the genotypes that
carry ML-accepted alleles, so an accepted allele is never emitted as `0/0`
because the prior-weighted posterior preferred reference; `0/0` is the
molecular genotype only when no allele passes. This is GATK's division of
labour between QUAL and the genotype likelihoods: the per-locus prior
(1e-3 for SNPs, 1e-4 for indels) and a library with frequent homopolymer
stutter would otherwise veto real heterozygotes under a passing ML score.
For molecular calls `GQ` is derived from the posterior error of the selected
genotype conditional on that choice, with an error floor at the maximum
posterior mass observed to be removed during pruning. This is confidence
conditional on the accepted alleles, the retained haplotypes, model parameters
and shared CpG rate estimates. It is not an externally calibrated error rate
or a bound on all omitted hypotheses.

Legacy `GL` and `GC` scores remain available on count-based calls and are omitted
for molecular genotypes. They are not substituted for a complete molecular
likelihood vector. On the count path, reference confidence uses ML threshold
margins and variant confidence uses read-count distributions.

For ML-filtered indels, VCF QUAL remains the forest's variant-presence confidence,
while GQ describes the genotype decision. A high QUAL therefore does not imply
high confidence in `1/1`. The two scores are also kept separate for compound
indel records.

## Strand-Specific Evidence

The count model avoids methylation-confounded observations when genotyping
`C→T` and `G→A` SNPs. Molecular genotyping uses the shared CpG shrinkage model
below, including reference and candidate de-novo CpGs, to account for conversion
in the read likelihood rather than treating it as independent genetic support.

## Phasing

`--phase` can be used to group heterozygous variants that share a read pair
into phase blocks and written with a `|` separator and a `PS` field:

```
chr12  20028875  .  A  G  ...  GT:PS  0|1:20028875
chr12  20028894  .  C  T  ...  GT:PS  1|0:20028875
```

Both records carry `PS=20028875`, the position of the block's first heterozygous site, so they belong to one block.
Within a block, the alleles written first all sit on the same chromosome copy: above, the reference `A` and the alternate `T` are on one copy, the alternate `G` and the reference `C` on the other.
Each block opens with its lower genotype allele first (`0|1` for a reference/alternate genotype, or e.g. `1|2` for two alternates); nothing relates one block to another, and heterozygous sites no read pair links stay unphased without `PS`.

Phasing now uses fragment likelihoods and a local PairHMM, conditioned on the two alleles of each **called** heterozygous SNP. This final phasing step changes chromosome assignments only; genotype selection occurs earlier using the model-specific procedure above. The previous signed-quality graph solver is used only as a reference in unit tests.

Full sequences, qualities, CIGAR coordinates and OT/OB orientation are captured before seqair recycles its read store. Read name plus read group identifies a molecule, and boundary alignments fetched again in byte-budget subdivisions are captured once. Overlapping mates form one consensus observation: agreement retains the better quality, while disagreement reduces it by the competing quality. Disagreeing CIGAR paths or reference skips use the aligned likelihood fallback.

A site needs at least two independent, quality-qualified molecules supporting **each** called allele before it can provide a genetic link. Consensus base and source mapping quality must be at least 20, and two molecules must link an edge. Converted-strand `T` at a C/T CpG and `A` at a G/A CpG cannot establish these links. Once a block is selected, its likelihood uses all eligible captured fragments covering its sites, including conversion-confounded observations and fragments reaching only one site. Adjacent candidate C/G partners are included together for CpG context even when they lack a genetic link.

Connected chromosome hypotheses expand progressively, with a beam of 32 phase states and 16 held-out states reconsidered on the next extension. Competing orientations at every selected site are protected; protection can enlarge the beam to 64 states. Each state pairs a chromosome with its complementary called alleles, so no alternate unphased genotypes are introduced. Retained states have equal phase priors. A nearby, independently supported CIGAR indel of at most 20 bp can be included as one nuisance allele, marginalizing its genotype and phase without emitting phased indels.

The neutral pilot and final PairHMM scores use the shared CpG shrinkage model described below: default conversion is 0.95, background conversion 0.003, with Beta(2,2) shrinkage and leave-one-fragment-out rate estimates. Homozygous SNP calls are applied to both chromosome backgrounds. Candidate sequence adjacency after SNP and indel edits determines reference and de-novo CpGs.

Each mate is aligned against its CIGAR-coordinate footprint padded by 12 bp plus the maximum total candidate indel length. Concordant overlaps use a joint molecular profile. An unsequenced gap between mates provides no alignment observations. Scores are reused for identical molecular profiles, template sequences, conversion rates, transition qualities and bands. One extra template base preserves CpG partners at window boundaries. Compair computes the forward sum with insertion emission one; no separate path-mass normalization is applied. Each molecule's relative evidence is capped by its mapping quality and the library profile's maximum (Q45 by default).

A genetic edge is emitted only if the retained-state posterior probability of its chosen relative phase is at least 0.99. Ambiguous edges can split a connected inference region into separate output phase sets. These probabilities are conditional on called genotypes, retained hypotheses and estimated shared rates, and are not externally calibrated phase quality scores. Genotype scores and the `GTM` method annotation survive this step; only `GT` orientation and `PS` are changed.

Use repeated `--phase-mask FILE` arguments to exclude centromeres and poorly mappable intervals from **phasing only**. Inputs can be plain, gzip or BGZF BED files, including the UCSC centromere table with its `Chr Start End GapType` header. Coordinates are zero-based and half-open; overlapping intervals from all files are merged. Masked positions retain variant and methylation calling but have no phased genotype or `PS`. For models using molecular features, masks also make those features unavailable and select count-based genotyping (`GTM=counts`); consequently calls can differ from an unmasked molecular run. Masks are barriers even when a read pair bridges them: inference, CpG context, nearby indel proposals and PairHMM flanks are restricted to each unmasked interval. Paths must be supplied explicitly, and contig names must match the BAM/reference.

```bash
rastair call --phase --phase-mask centromeres.txt --phase-mask blacklist.bed.gz \
  --fasta-file reference.fa reads.bam --vcf calls.vcf.gz
```

Indel records are not phased: the likelihood phaser links called heterozygous SNVs, and the pair-HMM indel caller's records (`--experimental-indels`, see [indels](indels.md)) carry no `PS`.

```admonish warning
Phasing cannot be combined with `--cpgs-only` or BED output alone. Inference is bounded to 32 heterozygous sites and a 2 kb span. Blocks also stop at processing-segment boundaries; overlap sites contribute evidence, but each phase set is anchored on its first emitted site. The production phasing model currently assumes TAPS chemistry with the defaults above.
```

`scripts/assess_phasing.py` compares called SNP phase with independent phased truth, such as the [Illumina Platinum Genomes NA12878 calls](https://github.com/illumina/platinumgenomes). Supply its confident-region BED and any masks with `--exclude-bed`. Genotype errors are reported separately; switch errors are orientation changes between adjacent genotype-matched phased SNPs sharing both query and truth phase sets. A chromosome-wide orientation flip is irrelevant, while an isolated erroneous SNP orientation contributes two switches. Singleton intersections cannot assess phase. Report phase coverage alongside error rates, since conservative thresholds can improve the rate by leaving difficult links unphased.

## The molecular likelihood model

The regional SNV features, molecular genotyping and phasing share one likelihood model over the captured read pairs. Every focal position is evaluated against the observed alleles at its site and its neighbours'; the regional search is described under [machine-learning models](ml.md#haplotype-likelihood-features). Conversion of methylated and unmethylated cytosines defaults to `0.95` and `0.003`; these are separate from base-quality sequencing errors, and the lambda spike-in is not used to fit them, since incomplete lambda conversion can underestimate conversion efficiency by 2–5 percentage points.

Both mates share a chromosome assignment, and overlapping bases form one molecular consensus observation and share a conversion event. Agreement retains the better base quality, without adding qualities. Disagreement retains the better-quality base and reduces its Phred quality by the competing quality: Q40 versus Q30 becomes Q10, while equal-quality conflicts become Q0 and carry no allele preference. Unknown bases carry no competing evidence. Mapping quality breaks ties between equal base qualities. A conversion-neutral pilot pass gives C/T identical emissions on OT at potential CpG C positions, and G/A identical emissions on OB at potential CpG G positions; potential positions are the union across candidate sequences, including de-novo and indel-created CpGs. The pilot supplies genotype and molecular chromosome probabilities for estimating one methylation rate per candidate CpG, shared across haplotypes that contain it. A symmetric Beta(2,2) prior shrinks low-information rates toward 0.5, only molecules with pilot CpG-carrier probability at least 0.9 can fit a rate, each fragment is scored with its own conversion observation's contribution removed, and candidate CpGs without usable observations retain the prior mean. These are approximate plug-in scores and do not model allele-specific methylation.

### Local PairHMM realignment

The forward algorithm has match, insertion and deletion states and runs on the banded SIMD kernels of [compair](https://github.com/Softleif/seqair/tree/main/crates/compair). It uses the GATK forward convention: insertion emission one, a uniform template-start prior, and no division by a separate read-length-dependent path mass. These are likelihood scores used to compare candidate sequences, not probabilities normalized over all possible read strings. Match emissions combine TAPS conversion and sequencing error; unknown template bases are uninformative.

Non-repeat gaps default to Phred 45 opening and Phred 10 extension. Repeat-dependent insertion opening, deletion opening and extension qualities come from the library's error profile, evaluated on the reference repeat context at the positions a read aligns to and projected into compair's read transition arrays using the original alignment coordinates; this is not a path-dependent transition model at every possible alignment column. Original coordinates set the band and context projection; the forward sum still includes alternative alignments.

Molecular profiles retain the sequenced bases within the unmasked inference domain. Each mate or joint overlap profile is aligned against its CIGAR-coordinate footprint padded by 12 bp plus the maximum total indel length of a candidate chromosome, in a band widened for observed and candidate indels; all compared chromosomes use the same reference interval, and soft clips are excluded. Paths outside the band remain an approximation. This is candidate-based regional realignment, without a de Bruijn graph assembler. For concordant overlapping CIGAR paths, mates form one molecular profile with a consensus observation at each shared base. Disjoint mates remain separate profiles assigned to the same chromosome. Discordant overlapping CIGAR paths and reference skips fall back to aligned SNP emissions at the observed bases, with the same consensus policy and the candidate's edited CpG context. A mate outside the unmasked reference domain contributes no sequence or orientation evidence.

The regional SNV model proposes no indels: a read carrying one is scored through its gap against chromosomes that differ by SNVs only. The phaser proposes the one nearest CIGAR insertion or deletion of at most 20 bp with two supporting reads as a nuisance allele, left-normalized against the local reference, with a genotype prior of `1 - 0.0001` for reference and `0.0001` spread over the alternate genotypes, and marginalizes over its presence and phase without emitting it. Candidate sequence adjacency determines CpGs after SNPs or that indel create or destroy them, and the shared shrinkage also estimates rates for CpGs in realignment flanks.

For each molecule, both the neutral CpG pilot and the final likelihood matrix apply the same symmetric floor: no candidate can score more than `min(MAPQ, maximum_evidence_phred)` below that molecule's best candidate, in Phred units. The cap applies once after combining the mates. Missing mapping quality or wholly uninformative observations neutralize the molecule. It limits outliers in either allele direction; it does not replace the repeat error model.
