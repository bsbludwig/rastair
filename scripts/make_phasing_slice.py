#!/usr/bin/env python3
"""Cut the NA12878 phasing-accuracy fixture out of the chr12 TAPS BAM.

Writes, into OUT_DIR, a 200 kb window of chr12 re-based onto its own contig:
`ref.fa(.fai)`, `reads.bam(.bai)`, `truth.vcf` (Platinum Genomes, phased) and
`confident.bed`. Only templates whose two primary reads both lie inside the
window are kept, so every mate link stays valid.

    scripts/make_phasing_slice.py tmp/taps tmp/na12878_phasing
    RASTAIR_PHASING_SLICE=tmp/na12878_phasing \
      cargo test --features experimental-seqair --test phasing_cli phasing_accuracy

The reads are partner data and stay out of git.

Needs samtools and bcftools on PATH, and in DATA_DIR: NA12878_aa_chr12.bam,
hg38.fa.gz, NA12878_PG_chr12.phased.vcf.gz and PG_ConfidentRegions_hg38.bed.gz.
"""
import re
import subprocess
import sys
from pathlib import Path

CHROM, START, END = "chr12", 31_800_000, 32_000_000  # 0-based, half-open
CONTIG = f"{CHROM}_{START}"
LENGTH = END - START


def run(*args, **kwargs):
    return subprocess.run(args, check=True, capture_output=True, text=True, **kwargs).stdout


def ref_span(cigar):
    return sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar) if op in "MDN=X")


def main(data, out):
    out.mkdir(parents=True, exist_ok=True)
    region = f"{CHROM}:{START + 1}-{END}"

    fasta = run("samtools", "faidx", str(data / "hg38.fa.gz"), region).splitlines()[1:]
    (out / "ref.fa").write_text(f">{CONTIG}\n" + "\n".join(fasta) + "\n")
    run("samtools", "faidx", str(out / "ref.fa"))

    sam = run("samtools", "view", "-F", "0x900", str(data / "NA12878_aa_chr12.bam"), region)
    records = [line.split("\t") for line in sam.splitlines()]
    inside = {}
    for fields in records:
        pos = int(fields[3]) - 1
        ok = fields[6] == "=" and pos >= START and pos + ref_span(fields[5]) <= END
        inside.setdefault(fields[0], []).append(ok)
    keep = {name for name, oks in inside.items() if len(oks) == 2 and all(oks)}

    header = run("samtools", "view", "-H", str(data / "NA12878_aa_chr12.bam")).splitlines()
    lines = [h for h in header if h.startswith(("@HD", "@RG"))]
    lines.insert(1, f"@SQ\tSN:{CONTIG}\tLN:{LENGTH}")
    for fields in records:
        if fields[0] not in keep:
            continue
        fields[2] = CONTIG
        fields[3] = str(int(fields[3]) - START)
        fields[7] = str(int(fields[7]) - START)
        # PG names a pipeline path, SA points at other contigs.
        fields = fields[:11] + [t for t in fields[11:] if not t.startswith(("PG:", "SA:"))]
        lines.append("\t".join(fields))
    run("samtools", "view", "-b", "-o", str(out / "reads.bam"), "-", input="\n".join(lines) + "\n")
    run("samtools", "index", str(out / "reads.bam"))

    confident = []
    for line in run("tabix", str(data / "PG_ConfidentRegions_hg38.bed.gz"), region).splitlines():
        chrom, a, b = line.split("\t")[:3]
        a, b = max(int(a), START) - START, min(int(b), END) - START
        if a < b:
            confident.append(f"{CONTIG}\t{a}\t{b}")
    (out / "confident.bed").write_text("\n".join(confident) + "\n")

    truth = run("bcftools", "view", "-H", "-v", "snps", "-r", region,
                str(data / "NA12878_PG_chr12.phased.vcf.gz"))
    body = []
    for line in truth.splitlines():
        fields = line.split("\t")
        fields[0] = CONTIG
        fields[1] = str(int(fields[1]) - START)
        fields[7] = "."
        body.append("\t".join(fields))
    (out / "truth.vcf").write_text(
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID={CONTIG},length={LENGTH}>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n"
        + "\n".join(body) + "\n"
    )
    print(f"{len(keep)} templates, {len(body)} truth SNVs, {len(confident)} confident intervals")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
