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

See [VCF Fields](vcf-fields.md) for a detailed description of the fields in the VCF output.
