# VCF

The main output of Rastair2 is a @VCF file.
It contains all metrics that Rastair2 calculates,
either for all @variant:pl, or only for @CpG sites.

# BCF output, compression

VCF files are text-based and can be quite large, especially for whole-genome sequencing data.
Rastair2 can also output @BCF files (binary VCF format) which are more compact and faster to read.
Alternatively, it can compress the VCF file transparently using @bgzip.
All formats can be read by `bcftools`, just like regular VCF files.

By specifying the file extension (`.vcf`, `.bcf`, or `.vcf.gz`) Rastair2 will automatically detect which format to write.
