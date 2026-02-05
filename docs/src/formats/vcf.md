# VCF

The main output of Rastair is a @VCF file.
It contains all metrics that Rastair calculates,
either for all @variant:pl, or only for @CpG sites.

## BCF output, compression

VCF files are text-based and can be quite large, especially for whole-genome sequencing data.
Rastair can also output @BCF files (binary VCF format) which are more compact and faster to read.
Alternatively, it can compress the VCF file transparently using @bgzip.
All formats can be read by `bcftools`, just like regular VCF files.

By specifying the file extension (`.vcf`, `.bcf`, or `.vcf.gz`) Rastair will automatically detect which format to write.

## Fields

By default, only a few fields are include.
You can enable more using the `--vcf-info-fields` and `--vcf-info-fields` flags with comma-separated field names.
See [VCF Fields](vcf-fields.md) for a detailed description of the fields in the VCF output.

### Methylation information
The estimated $ \beta $ value is encoded in the VCF file through the `M5mC` FORMAT field which is provided at all positions where `CPG` or `CPGnovo` are set. We include all CpG positions in the reference, whether they overlap a variant or not. Where they do not overlap a variant, the `ALT` is set to `.`, _aka_ the missing value.

While rare, it is possible to have more than one methylation call at one position. This can occur in `CCG` positions where the middle C is also a C>G SNP. This leads to two alleles, one where the middle position is a G and the other where it is a C. Each could be methylated at different levels. In these cases, the `M5mC` value will be a comma-separated list.

### Quality scores
Quality scores in the vcf output reflect the ML-derived probability that a variant is "true". When we train the ML model, we calibrate it against a validation dataset and estimate the fraction of true SNPs called at every score cutoff.
Based on this, we can fit a regression (["Platt scaling"](https://en.wikipedia.org/wiki/Platt_scaling)) that converts the score into a probability $ p $ that a SNP call is true. For the quality score, we encode the inverse as Phred: $ -10*log_{10}(1-p) $. For methylated positions where no variants were observed, we record the maximum Phred, _ie._ 99.