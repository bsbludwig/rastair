# Indel Calling

As of version 2.2.0, rastair supports calling @indel:pl variants in the `call` workflow.

```admonish warning
This is still an experimental feature!
The `--experimental-indels` command-line flag must be used to enable it.
```

## How Indel Candidates Are Built

At each pileup position, Rastair collects insertion and deletion observations from aligned reads, one per read pair, after the usual base-quality, mapping-quality, flag and overlap filters.
Each observed indel is represented as an allele, and the reads whose alignment carries it are counted.
These counts only decide which alleles are looked at; everything a record says about them comes from the realignment described next.
See the [CLI reference][cli-indels] for the options.

[cli-indels]: ../cli.md#indel-options

## How Indels Are Called

Indel positions at most 20 bases apart are called together, as one cluster: an aligner often spells one event differently in different reads, and one chromosome can carry two events close together.
A cluster's candidates are its indel alleles with at least `--min-indel-ao` reads, the best supported first, at most three per cluster and at most two at one position.
Rastair builds a local haplotype for every combination of candidates one chromosome can carry: the reference, each candidate, and candidates at different positions together (two at the same position, or a deletion over another candidate's position, cannot share a chromosome).
Every read covering the cluster is then aligned to each haplotype by a pair hidden Markov model that knows about TAPS conversion, so a converted base costs no more than it should.
Gaps are cheaper inside the read's own tandem repeats, where polymerase slippage makes them common.
Only one read of each read pair is used (the one reaching further), and only reads with a mapping quality of at least the `call` minimum and a known strand.

Positions near the cluster where at least five reads show a base other than the reference, one that TAPS conversion cannot explain, are left open in every haplotype and match any base.
Otherwise a nearby SNV, missing from every haplotype, would count against the reads that carry it, and an indel that happens to line those reads up would be called in its place.

From these alignments, Rastair finds the pair of haplotypes, one per chromosome copy, that best explains the reads.
A candidate is called when the reads are at least ten times likelier with it present than without it, comparing the best pair carrying it to the best account that does not: a pair without the candidate, or such a pair with a small share of its reads (10 %) drawn from the candidate as an artefact, the way stutter, a mis-mapped paralog or a contaminant would supply them, that account paying a hundredfold prior.
So seven reads of an allele among seventy-three are an artefact, not a genotype, while four among eight are a het.
That ratio, as a Phred score, is the record's `QUAL`, so a called indel always has a `QUAL` of at least 10.
The genotype is the likeliest pair carrying exactly the called candidates, normally the best pair itself, and each record's genotype counts the chromosomes of that pair carrying its allele: two deletions 13 bases apart that every read carries are each `1/1`.
`GC` is the confidence that the written genotype is right, over all pairs, so it is low whenever the written pair is not the best one.
Two called alleles at one position make one record with a `1/2` genotype.
`DP` is the number of reads the realignment assigns to a haplotype, and `AD` how many of them it assigns to a haplotype without the record's allele and with it.
A read is assigned to its likeliest haplotype when that is at least 1.6 times (0.2 log10) likelier than the runner-up: a lower bar than calling, because in a long repeat one read says only a little about where the gap sits.
A read that fits two haplotypes about equally, as one ending before the indel does, counts towards neither `AD` nor `DP`.
On chr12 of NA12878 (~26x) against Platinum Genomes, this reached an F1 of 0.905 on insertions and 0.923 on deletions; at ~73x against GIAB, 0.971 and 0.980.

There is no strand filter: an allele seen on one bisulfite strand only is called like any other.
Under `--all`, an allele the pair-HMM looked at but did not call is written as `0/0` with the `low_indel_qual` filter.

Indel calling needs a build with the `experimental-seqair` backend, since it aligns reads while their position is built; other builds refuse `--experimental-indels`.
No @ML model is involved, so `--ml` and `--no-ml` do not affect indel calls.
With `--phase`, heterozygous indels are phased along with the SNVs (see [genotyping](genotyping.md#phasing)).
