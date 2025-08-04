# De-novo CpGs

In addition to looking at methylation of known @CpG sites,
i.e. sites that are already present in the @referenceGenome,
Rastair can also call @denovo sites.
These are sites where a `C` and/or `G` @altAllele exists in @read:pl
followed/preceded by a `G` and/or `C` @refAllele,
and thus creating new CpG sites.

## Methylation of de-novo CpGs

After Rastair has identified de-novo CpG sites,
it will also call methylation for these sites.
For a @methylated de-novo CpG site,
there have to be both `C` and `T` (or `G` and `A`) @altAllele:pl,
which means the amount of evidence present is generally lower than for known CpG sites.

The same filter criteria as for known CpG sites is applied.
