#!/usr/bin/env python3
"""Compare haplotype likelihood engines on matched TAPS and ordinary WGS libraries.

Input JSON keys: reference, truth, confidence_bed, gatk, taps, wgs, region.
Truth is used after inference, only within confident regions. Candidate-based
calibration statistics are not whole-genome accuracy measurements.
"""

import argparse
import csv
import json
import math
import pathlib
import re
import subprocess


def query(tabix, path, region):
    return subprocess.check_output([tabix, str(path), region], text=True).splitlines()


def variants(lines):
    snps, indels = set(), []
    for line in lines:
        fields = line.split("\t")
        if fields[6] not in ("PASS", "."):
            continue
        reference, alternates = fields[3], fields[4].split(",")
        genotype = fields[9].split(":")[fields[8].split(":").index("GT")]
        indices = {int(a) for a in re.split(r"[/|]", genotype) if a not in (".", "")}
        for index in indices - {0}:
            alternate, position = alternates[index - 1], int(fields[1])
            if len(reference) == len(alternate) == 1 and alternate in "ACGT":
                snps.add((fields[0], position, reference, alternate))
            else:
                indels.append((fields[0], position, position + len(reference) - 1))
    return snps, indels


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))


def key(row):
    return row["chrom"], int(row["pos"]), row["ref"], row["alt"]


def score(data, truth, callable_site):
    data = [row for row in data if callable_site(row["chrom"], int(row["pos"]))]
    evaluated = {key(row) for row in data}
    calls = {key(row) for row in data if float(row["posterior"]) >= 0.5}
    tp, fp, fn = len(calls & truth), len(calls - truth), len(truth - calls)
    loss, brier = 0.0, 0.0
    for row in data:
        p = min(max(float(row["posterior"]), 1e-15), 1.0 - 1e-15)
        label = key(row) in truth
        loss -= math.log(p if label else 1.0 - p)
        brier += (p - int(label)) ** 2
    return {
        "evaluated_alleles": len(data), "truth_snps": len(truth),
        "truth_snps_with_evidence": len(evaluated & truth),
        "tp": tp, "fp": fp, "fn": fn,
        "precision": tp / (tp + fp) if tp + fp else None,
        "recall": tp / (tp + fn) if tp + fn else None,
        "mean_log_loss_on_candidates": loss / len(data) if data else None,
        "brier_score_on_candidates": brier / len(data) if data else None,
        "cpg_context_limited_rows": sum(row["cpg_context_limited"] == "true" for row in data),
        "alignment_indel_rows": sum(int(row["alignment_indel_reads"]) > 0 for row in data),
        "realigned_rows": sum(int(row.get("realigned_fragments", 0)) > 0 for row in data),
        "aligned_fallback_rows": sum(int(row.get("aligned_fallback_fragments", 0)) > 0 for row in data),
        "omitted_indel_rows": sum(int(row.get("omitted_indels", 0)) > 0 for row in data),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=pathlib.Path, required=True)
    parser.add_argument("--rastair", type=pathlib.Path, required=True)
    parser.add_argument("--tabix", default="tabix")
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    parser.add_argument("--indel-padding", type=int, default=10)
    parser.add_argument("--pair-hmm", action="store_true", help="also evaluate TAPS-aware and ordinary local PairHMM arms")
    parser.add_argument("--include-indel-neighborhoods", action="store_true", help="score SNPs near truth indels and report them separately")
    args = parser.parse_args()
    if args.indel_padding < 0:
        parser.error("indel padding must be nonnegative")
    inputs = json.loads(args.inputs.read_text())
    for field in ["reference", "truth", "confidence_bed", "gatk", "taps", "wgs"]:
        if not pathlib.Path(inputs[field]).is_file():
            parser.error(f"missing {field}: {inputs[field]}")
    region = inputs["region"]
    truth, indels = variants(query(args.tabix, inputs["truth"], region))
    gatk, _ = variants(query(args.tabix, inputs["gatk"], region))
    confident = [line.split("\t") for line in query(args.tabix, inputs["confidence_bed"], region)]

    def near_indel(chrom, pos):
        return any(
            c == chrom and start - args.indel_padding <= pos <= end + args.indel_padding for c, start, end in indels
        )

    def callable_site(chrom, pos):
        return any(c == chrom and int(start) <= pos - 1 < int(end) for c, start, end, *_ in confident) and (
            args.include_indel_neighborhoods or not near_indel(chrom, pos)
        )

    truth = {allele for allele in truth if callable_site(allele[0], allele[1])}
    gatk = {allele for allele in gatk if callable_site(allele[0], allele[1])}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs, commands = {}, {}
    discovered = set()
    for library, non_taps in [("wgs", True), ("taps", False)]:
        output = args.output_dir / f"{library}_discovery.tsv"
        command = [str(args.rastair.resolve()), "ml", "haplotypes", "--fasta-file", inputs["reference"], inputs[library],
                   "--region", region, "--min-alt-observations", "2", "--max-sites", "1", "--output", str(output)]
        if non_taps:
            command.append("--non-taps")
        with (args.output_dir / f"{library}_discovery.log").open("w") as log:
            subprocess.run(command, check=True, stdout=log, stderr=log)
        discovered.update((row["chrom"], int(row["pos"])) for row in rows(output))
        commands[f"{library}_discovery"] = command
    discovered.update((chrom, pos) for chrom, pos, *_ in truth)
    targets = args.output_dir / "sites.bed"
    targets.write_text("".join(f"{chrom}\t{pos - 1}\t{pos}\n" for chrom, pos in sorted(discovered)))
    arms = [
        ("wgs", "wgs", True, 16, "aligned"), ("taps", "taps", False, 16, "aligned"),
        ("taps_non_taps_model", "taps", True, 16, "aligned"), ("taps_quadrature_32", "taps", False, 32, "aligned"),
    ]
    if args.pair_hmm:
        arms.extend([
            ("wgs_pair_hmm", "wgs", True, 16, "pair-hmm"), ("taps_pair_hmm", "taps", False, 16, "pair-hmm"),
            ("taps_pair_hmm_non_taps_model", "taps", True, 16, "pair-hmm"),
            ("taps_pair_hmm_quadrature_32", "taps", False, 32, "pair-hmm"),
        ])
    for name, library, non_taps, order, engine in arms:
        output = args.output_dir / f"{name}.tsv"
        command = [str(args.rastair.resolve()), "ml", "haplotypes", "--fasta-file", inputs["reference"], inputs[library],
                   "--region", region, "--conversion-rate", "0.95", "--false-conversion-rate", "0.003",
                   "--min-alt-observations", "2", "--quadrature-order", str(order), "--likelihood-engine", engine,
                   "--sites-bed", str(targets), "--output", str(output)]
        if non_taps:
            command.append("--non-taps")
        with (args.output_dir / f"{name}.log").open("w") as log:
            subprocess.run(command, check=True, stdout=log, stderr=log)
        commands[name], outputs[name] = command, rows(output)
    metrics = {name: score(data, truth, callable_site) for name, data in outputs.items()}
    metrics["gatk_wgs"] = {"tp": len(gatk & truth), "fp": len(gatk - truth), "fn": len(truth - gatk)}
    near_indel_metrics = {
        name: score(data, {allele for allele in truth if near_indel(allele[0], allele[1])},
                    lambda chrom, pos: callable_site(chrom, pos) and near_indel(chrom, pos))
        for name, data in outputs.items()
    } if args.include_indel_neighborhoods else None
    coarse = {key(row): float(row["posterior"]) for row in outputs["taps"]}
    fine = {key(row): float(row["posterior"]) for row in outputs["taps_quadrature_32"]}
    if coarse.keys() != fine.keys():
        raise ValueError("quadrature comparison changed the candidate set")
    difference = max((abs(coarse[k] - fine[k]) for k in coarse), default=0.0)
    hmm_difference = None
    if args.pair_hmm:
        hmm_coarse = {key(row): float(row["posterior"]) for row in outputs["taps_pair_hmm"]}
        hmm_fine = {key(row): float(row["posterior"]) for row in outputs["taps_pair_hmm_quadrature_32"]}
        if hmm_coarse.keys() != hmm_fine.keys():
            raise ValueError("PairHMM quadrature comparison changed the candidate set")
        hmm_difference = max((abs(hmm_coarse[k] - hmm_fine[k]) for k in hmm_coarse), default=0.0)
    report = {
        "models": ["aligned_snp_v2"] + (["local_pair_hmm_v2"] if args.pair_hmm else []), "inputs": inputs, "commands": commands,
        "scoring": "SNP allele presence at common focal positions in GIAB confident regions",
        "include_indel_neighborhoods": args.include_indel_neighborhoods,
        "indel_padding": args.indel_padding, "metrics": metrics,
        "near_indel_metrics": near_indel_metrics,
        "quadrature_max_posterior_difference": difference,
        "pair_hmm_quadrature_max_posterior_difference": hmm_difference,
        "lambda_policy": "Diagnostic only; partial lambda conversion can underestimate efficiency by 2-5 percentage points",
    }
    (args.output_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"metrics": metrics, "quadrature_max_posterior_difference": difference,
                      "pair_hmm_quadrature_max_posterior_difference": hmm_difference}, indent=2))


if __name__ == "__main__":
    main()
