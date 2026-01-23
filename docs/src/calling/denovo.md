# De-novo CpGs

In addition to looking at methylation of known @CpG sites,
i.e. sites that are already present in the @referenceGenome,
Rastair can also call @denovo sites.
These are sites where a `C` and/or `G` @altAllele exists in @read:pl
followed/preceded by a `G` and/or `C` @refAllele,
and thus creating new CpG sites.

## Example

Here is an example of a de-novo CpG site in a pileup:

```
Position    1 2 3 4 5 6 7 8

Reference:  A T C C T A G C  Strand

Reads:      A T T G T A G C    +
            A T C A T A G C    -
            A T C A T A G C    -
            A T T G T A G C    +
            A T T G T A G C    +
                  ↑
                  De-novo CpG created by C>G variant
```

In this example, some reads have a `G` at position 3 where the reference has a `C`.
This creates a new CpG dinucleotide (CG) that is not present in the reference genome.
This newly generated CpG is methylated, which means that both the C and the G position
will show T/A at @OT/@OB reads, respectively!

## Methylation of de-novo CpGs

After Rastair has identified de-novo CpG sites,
it will also call methylation for these sites.
For a @methylated de-novo CpG site,
there have to be both `C` and `T` (or `G` and `A`) @altAllele:pl,
which means the amount of evidence present is generally lower than for known CpG sites.

The same filter criteria as for known CpG sites are applied.