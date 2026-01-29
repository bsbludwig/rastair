# Machine Learning Models

Rastair uses @ML to distinguish true @methylation calls and variants from sequencing artifacts (and noise).
The ML models evaluate each alt at each position and assign a prediction score for that alt to be a true variant.

By default, ML is enabled with a threshold of `0.8`.
A pre-trained model is bundled with Rastair.

The model handles three contexts:
- **@CpG methylation sites**: Standard 5mC detection in @CpG sites
- **@Denovo**: New methylation sites not in reference
- **Other variants**: Non-CpG @SNP:pl and @indel:pl

Some features are shared across models, but some are model specific. Here is a plot of all features used in each model, ranked by importance:

![Feature importance ranked by mean importance across models](../img/feature_importance.png "Feature importance per model")

Here, `_adj` refers to a feature of the `adjacent` nucleotide, _ie_ the C when evaluating a G or the G when evaluating a C. All scores that refer to allele counts are normalised either implicitly (where the score itself is a ratio) or explicitly to the total depth at that position. Of note, the `alt_score` is a simple ratiometric score to establish the base-quality weighted enrichment of variant reads over non-variant reads: $$ altscore = \frac{log_{2}(sb_{alt}*bq_{alt}+1)}{log_{2}(sb_{ref}*bq_{ref}+1)} $$, where $sb_{ref/alt}$ refers to the number of reads with ref/alt on the OT (for C positions) or OB (for G positions), and $bq_{ref/alt}$ is the corresponding strand-specific @RMS of base qualities.

## Adjusting the Threshold

Default:

```bash
rastair call input.bam
```

Change threshold:

```bash
rastair call --ml 0.9 input.bam
```

Disable ML (faster, less accurate):

```bash
rastair call --no-ml input.bam
```

## Performance Considerations

ML inference is the slowest part of the pipeline.
To speed up calling:

- **Reduce positions evaluated**: Use `--cpg-only` if you only care about CpG methylation
- **Disable ML**: Use `--no-ml` for quick exploratory runs where accuracy isn't critical

## Training Custom Models

```admonish warning
The `ml train` subcommand is insufficiently tested at this point, and might not produce optimal models!
```

If you have ground-truth data (e.g., validated @SNP:pl), you can train a custom model tailored to your specific dataset, coverage, or sample type.

Basic training command:

```bash
rastair ml train input.bam --truth ground_truth.vcf --output ./my_models
```

See [the CLI docs](../cli.md#rastair-ml-train)
or run `rastair ml train --help`
for an overview of all training options.
