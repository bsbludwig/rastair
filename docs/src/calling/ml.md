# Machine Learning Models

Rastair uses @ML to distinguish true @methylation calls and variants from sequencing artifacts (and noise).
The ML models evaluate each alt at each position and assign a prediction score for that alt to be a true variant.

By default, ML is enabled with a threshold of `0.8`.
A pre-trained model is bundled with Rastair.

The model handles three contexts:
- **@CpG methylation sites**: Standard 5mC detection in @CpG sites
- **@Denovo**: New methylation sites not in reference
- **Other variants**: Non-CpG @SNP:pl and @indel:pl

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

If you have ground-truth data (e.g., validated @SNP:pl), you can train a custom model tailored to your specific dataset, coverage, or sample type.

Basic training command:

```bash
rastair ml train input.bam --truth ground_truth.vcf --output ./my_models
```

See [the CLI docs](../cli.md#rastair-ml-train)
or run `rastair ml train --help`
for an overview of all training options.
