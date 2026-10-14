#!/usr/bin/env python3
"""Assess candidate discovery, calibration and context loss on predefined windows.

The manifest contains inputs (as in validate_haplotype_likelihoods.py), regions
with region/split fields, negative_sites_per_region, seed and optional known_sites.
All truth SNP positions and near-indel candidate positions are retained. Other
negative candidate positions are sampled deterministically, with inverse sampling
weights for calibration. Accuracy counts describe the evaluated panel, not all
negative positions in the windows. No caller parameter is fitted here.
"""

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import pathlib
import re
import shutil
import subprocess
import time

from validate_haplotype_likelihoods import key, query, rows, variants

GENOTYPES = [(a, b) for i, b in enumerate("ACGT") for a in "ACGT"[:i + 1]]
ARMS = [("wgs_aligned", "wgs", "aligned"), ("taps_aligned", "taps", "aligned"),
        ("wgs_pair_hmm", "wgs", "pair-hmm"), ("taps_pair_hmm", "taps", "pair-hmm")]


def arms(manifest):
    configured = manifest.get("arms")
    if configured is None:
        return ARMS
    result = [(a["name"], a["library"], a["engine"]) for a in configured]
    if not result or len({a[0] for a in result}) != len(result) or any(
            not re.fullmatch(r"[a-z0-9_]+", name) or library not in ("wgs", "taps") or
            engine not in ("aligned", "pair-hmm") for name, library, engine in result):
        raise ValueError("Invalid benchmark arms")
    return result


def arm_arguments(manifest, name):
    return next((a.get("arguments", []) for a in manifest.get("arms", []) if a["name"] == name), [])


def digest(path):
    result = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def truth_genotypes(lines):
    result = {}
    for line in lines:
        f = line.split("\t")
        if f[6] not in ("PASS", "."):
            continue
        if len(f[3]) != 1 or f[3] not in "ACGT":
            continue
        bases = [f[3], *f[4].split(",")]
        gt = f[9].split(":")[f[8].split(":").index("GT")]
        alleles = re.split(r"[/|]", gt)
        if len(alleles) != 2 or any(a == "." for a in alleles):
            continue
        called = [bases[int(a)] for a in alleles]
        if all(a in "ACGT" and len(a) == 1 for a in called):
            site = f[0], int(f[1])
            genotype = tuple(sorted(called))
            if site in result and result[site] != genotype:
                raise ValueError(f"Conflicting truth genotypes at {site}")
            result[site] = genotype
    return result


def select_panel(candidates, truth_sites, mandatory, count, seed):
    pool = sorted(candidates - truth_sites - mandatory)
    sampled = sorted(pool, key=lambda k: hashlib.sha256(f"{seed}:{k[0]}:{k[1]}".encode()).digest())[:count]
    weight = len(pool) / len(sampled) if sampled else 1.0
    result = {site: {"weight": 1.0, "stratum": "truth"} for site in truth_sites}
    result.update({site: {"weight": 1.0, "stratum": "mandatory_negative"} for site in mandatory - truth_sites})
    result.update({site: {"weight": weight, "stratum": "sampled_negative"} for site in sampled})
    return result, {"negative_pool": len(pool), "negative_sample": len(sampled), "negative_weight": weight}


def run(command, prefix):
    resource_path = prefix.with_suffix(".resources.tsv")
    started = time.monotonic()
    with prefix.with_suffix(".log").open("w") as log:
        subprocess.run(["/usr/bin/time", "-f", "%e\t%M", "-o", str(resource_path), *command],
                       check=True, stdout=log, stderr=log)
    elapsed, rss = resource_path.read_text().strip().split("\t")
    return {"command": command, "wall_seconds": time.monotonic() - started,
            "tool_seconds": float(elapsed), "peak_rss_kib": int(rss)}


def region_data(inputs, region, tabix):
    lines = query(tabix, inputs["truth"], region)
    truth, indels = variants(lines)
    intervals = [(c, int(a), int(b)) for c, a, b, *_ in
                 (line.split("\t") for line in query(tabix, inputs["confidence_bed"], region))]
    def confident(site):
        return any(c == site[0] and a <= site[1] - 1 < b for c, a, b in intervals)
    truth = {k for k in truth if confident(k[:2])}
    gt = {k: v for k, v in truth_genotypes(lines).items() if confident(k)}
    return truth, gt, indels, intervals, confident


def prepare_region(selected, manifest, root, binary, tabix):
    region = selected["region"]
    chrom = region.split(":")[0]
    folder = root / chrom
    folder.mkdir(exist_ok=True)
    inputs = manifest["inputs"]
    truth, gt, indels, intervals, confident = region_data(inputs, region, tabix)
    candidates, commands, discovery = set(), {}, {}
    for name, library, engine in arms(manifest):
        output = folder / f"{name}_discovery.tsv"
        command = [str(binary), "ml", "haplotypes", "--fasta-file", inputs["reference"], inputs[library],
                   "--region", region, "--min-alt-observations", "2", "--likelihood-engine", engine,
                   "--candidates-only", "--output", str(output)]
        if library == "wgs":
            command.append("--non-taps")
        commands[name] = run(command, folder / f"{name}_discovery")
        sites = {(r["chrom"], int(r["pos"])) for r in rows(output) if confident((r["chrom"], int(r["pos"])))}
        truth_sites = {k[:2] for k in truth}
        discovery[name] = {"confident_candidates": len(sites), "truth_sites": len(truth_sites),
                           "truth_sites_discovered": len(sites & truth_sites),
                           "missed_truth_sites": sorted(truth_sites - sites)}
        candidates.update(sites)
    truth_sites = {k[:2] for k in truth}
    def near(site):
        return any(c == site[0] and a - 10 <= site[1] <= b + 10 for c, a, b in indels)
    mandatory = {site for site in candidates if near(site)}
    mandatory.update(tuple(site) for site in manifest.get("known_sites", []) if site[0] == chrom and confident(site))
    panel, sampling = select_panel(candidates, truth_sites, mandatory,
                                   manifest["negative_sites_per_region"], manifest["seed"])
    bed = folder / "sites.bed"
    bed.write_text("".join(f"{c}\t{p - 1}\t{p}\n" for c, p in sorted(panel)))
    match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", region)
    if not match:
        raise ValueError(f"Invalid region: {region}")
    start, end = int(match[2]), int(match[3])
    prepared = {"selection": selected, "sampling": sampling, "panel_sites": len(panel),
                "panel": [{"chrom": c, "pos": p, **meta} for (c, p), meta in sorted(panel.items())],
                "discovery": discovery, "discovery_commands": commands,
                "candidate_union_sites": len(candidates), "truth_snp_alleles": len(truth),
                "forced_known_sites_outside_discovery": sorted(mandatory - candidates),
                "truth_snp_sites": len(truth_sites), "truth_indel_records": len(indels),
                "truth_snp_alleles_near_indels": sum(near(k[:2]) for k in truth),
                "bases": end - start + 1,
                "confident_bases": sum(max(0, min(end, b) - max(start - 1, a)) for c, a, b in intervals),
                "sites_sha256": digest(bed)}
    (folder / "prepared.json").write_text(json.dumps(prepared, indent=2) + "\n")
    print(f"Prepared {region}: {len(truth)} truth alleles, {len(candidates)} candidate positions, {len(panel)} panel positions", flush=True)
    return prepared


def reference_sequence(path, region):
    chrom, coordinates = region.split(":")
    lo, hi = map(int, coordinates.split("-"))
    index = {f[0]: tuple(map(int, f[1:5])) for f in
             (line.split("\t") for line in pathlib.Path(str(path) + ".fai").read_text().splitlines())}
    length, offset, width, line_bytes = index[chrom]
    start, end = max(0, lo - 2), min(length, hi + 1)
    with pathlib.Path(path).open("rb") as handle:
        handle.seek(offset + start // width * line_bytes + start % width)
        remaining = end - start
        sequence = bytearray()
        while remaining:
            line = handle.readline().strip()
            if not line:
                raise ValueError("Unexpected end of FASTA")
            chunk = line[:remaining]
            sequence.extend(chunk)
            remaining -= len(chunk)
    return start, sequence.decode().upper()


def validate(data, expected_sites):
    seen, sites = set(), {}
    for row in data:
        k = key(row)
        if k in seen:
            raise ValueError(f"Duplicate row {k}")
        seen.add(k)
        sites.setdefault(k[:2], []).append(row)
        for field in ["posterior", "p_het", "p_hom_alt", "expected_fraction", "site_posterior"]:
            value = float(row[field])
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"Invalid {field} at {k}")
        if abs(float(row["posterior"]) - float(row["p_het"]) - float(row["p_hom_alt"])) > 5e-12:
            raise ValueError(f"Invalid allele projection at {k}")
        if abs(float(row["expected_fraction"]) - float(row["p_het"]) / 2 - float(row["p_hom_alt"])) > 5e-12:
            raise ValueError(f"Invalid allele fraction at {k}")
        for field in ["genotype_posteriors", "site_genotype_posteriors"]:
            probabilities = list(map(float, row[field].split(",")))
            if len(probabilities) != 10 or any(not math.isfinite(p) or not 0 <= p <= 1 for p in probabilities) or abs(sum(probabilities) - 1) > 5e-12:
                raise ValueError(f"Invalid genotype distribution at {k}")
            present = sum(p for gt, p in zip(GENOTYPES, probabilities) if row["alt"] in gt)
            target = "posterior" if field == "genotype_posteriors" else "site_posterior"
            if abs(present - float(row[target])) > 5e-12:
                raise ValueError(f"Inconsistent allele/genotype posteriors at {k}")
        if int(row["retained_linking_fragments"]) > int(row["offered_linking_fragments"]):
            raise ValueError(f"Invalid linkage counts at {k}")
        if "read_observations" in row:
            observations = int(row["observations"])
            raw = int(row["read_observations"])
            overlap = int(row["overlapping_fragments"])
            discordant = int(row["discordant_fragments"])
            if not (0 <= discordant <= overlap <= observations <= int(row["fragments"])) or raw != observations + overlap:
                raise ValueError(f"Invalid molecular consensus counts at {k}")
    if set(sites) != expected_sites or any(len(v) != 3 for v in sites.values()):
        raise ValueError("Missing or extra target alleles")
    return len(data)


def decorate(data, prepared, truth, gt, indels, reference):
    panel = {(r["chrom"], r["pos"]): r for r in prepared["panel"]}
    start, sequence = reference
    for row in data:
        k = key(row)
        position = k[1] - 1 - start
        before = sequence[position - 1] if position else "N"
        after = sequence[position + 1] if position + 1 < len(sequence) else "N"
        ref_cpg = (k[2] == "C" and after == "G") or (k[2] == "G" and before == "C")
        denovo = (k[3] == "C" and after == "G") or (k[3] == "G" and before == "C")
        row.update(weight=panel[k[:2]]["weight"], label=k in truth,
                   context="reference_cpg" if ref_cpg else "single_snp_created_cpg" if denovo else "ordinary",
                   near_indel=any(c == k[0] and a - 10 <= k[1] <= b + 10 for c, a, b in indels),
                   truth_genotype=gt.get(k[:2]), split=prepared["selection"]["split"])
    return data


def allele_metrics(data, field="posterior"):
    tp = fp = fn = 0
    weight = loss = brier = weighted_fp = 0.0
    bins = [{"lo": a, "hi": b, "rows": 0, "weight": 0.0, "p_sum": 0.0, "truth_sum": 0.0}
            for a, b in [(0, .001), (.001, .01), (.01, .1), (.1, .5), (.5, .9), (.9, .99), (.99, 1.000000000001)]]
    for row in data:
        p, y, w = float(row[field]), int(row["label"]), row["weight"]
        called = p >= .5
        tp += int(called and y)
        fp += int(called and not y)
        fn += int(not called and y)
        weighted_fp += w * int(called and not y)
        clipped = min(max(p, 1e-15), 1 - 1e-15)
        weight += w
        loss -= w * math.log(clipped if y else 1 - clipped)
        brier += w * (p - y) ** 2
        bucket = next(b for b in bins if b["lo"] <= p < b["hi"])
        bucket["rows"] += 1
        bucket["weight"] += w
        bucket["p_sum"] += w * p
        bucket["truth_sum"] += w * y
    for b in bins:
        b["mean_probability"] = b["p_sum"] / b["weight"] if b["weight"] else None
        b["observed_fraction"] = b["truth_sum"] / b["weight"] if b["weight"] else None
    return {"evaluated_alleles": len(data), "tp": tp, "fp": fp, "fn": fn,
            "panel_precision": tp / (tp + fp) if tp + fp else None,
            "recall": tp / (tp + fn) if tp + fn else None,
            "estimated_candidate_population_fp": weighted_fp,
            "weighted_alleles": weight, "weighted_log_loss": loss / weight if weight else None,
            "weighted_brier": brier / weight if weight else None, "reliability_bins": bins}


def genotype_metrics(data, field="genotype_posteriors"):
    sites = {key(row)[:2]: row for row in data if row["truth_genotype"] is not None}
    errors, loss, brier = [], 0.0, 0.0
    for site, row in sites.items():
        probabilities = list(map(float, row[field].split(",")))
        target = GENOTYPES.index(row["truth_genotype"])
        best = max(range(10), key=lambda i: probabilities[i])
        loss -= math.log(max(probabilities[target], 1e-15))
        brier += sum((p - int(i == target)) ** 2 for i, p in enumerate(probabilities))
        if best != target:
            errors.append({"chrom": site[0], "pos": site[1], "truth": row["truth_genotype"],
                           "called": GENOTYPES[best], "posterior": probabilities[best]})
    return {"truth_snp_sites": len(sites), "genotype_errors": len(errors),
            "accuracy": 1 - len(errors) / len(sites) if sites else None,
            "mean_log_loss": loss / len(sites) if sites else None,
            "mean_brier": brier / len(sites) if sites else None, "errors": errors}


def summarize(data):
    sites = list({key(r)[:2]: r for r in data}.values())
    limited = [r for r in sites if r.get("cpg_context_limited", "false") == "true"]
    truth_sites = [r for r in sites if r["truth_genotype"] is not None]
    offered = sum(int(r["offered_linking_fragments"]) for r in sites)
    retained = sum(int(r["retained_linking_fragments"]) for r in sites)
    site_weight = sum(r["weight"] for r in sites)
    offered_weight = sum(r["weight"] * int(r["offered_linking_fragments"]) for r in sites)
    lost_weight = sum(r["weight"] * (int(r["offered_linking_fragments"]) - int(r["retained_linking_fragments"])) for r in sites)
    return {"alleles": allele_metrics(data), "genotypes": genotype_metrics(data),
            "site_only_alleles": allele_metrics(data, "site_posterior"),
            "site_only_genotypes": genotype_metrics(data, "site_genotype_posteriors"),
            "strata": {name: allele_metrics([r for r in data if predicate(r)]) for name, predicate in [
                ("reference_cpg", lambda r: r["context"] == "reference_cpg"),
                ("single_snp_created_cpg", lambda r: r["context"] == "single_snp_created_cpg"),
                ("ordinary", lambda r: r["context"] == "ordinary"),
                ("near_indel", lambda r: r["near_indel"]),
                ("limited_context", lambda r: r.get("cpg_context_limited", "false") == "true"),
                ("depth_under_10", lambda r: int(r["fragments"]) < 10),
                ("depth_10_to_19", lambda r: 10 <= int(r["fragments"]) < 20),
                ("depth_20_or_more", lambda r: int(r["fragments"]) >= 20)]},
            "context": {"sites": len(sites), "limited_sites": len(limited),
                "disconnected_sites": sum(int(r.get("disconnected_sites", 0)) for r in sites),
                "resource_limited_sites": sum(int(r.get("resource_limited_sites", 0)) for r in sites),
                "pruned_haplotypes": sum(int(r.get("pruned_haplotypes", 0)) for r in sites),
                "rescued_haplotypes": sum(int(r.get("rescued_haplotypes", 0)) for r in sites),
                "max_context_sites": max((int(r.get("context_sites", 1)) for r in sites), default=0),
                "maximum_pruned_mass": max((float(r.get("maximum_pruned_mass", 0)) for r in sites), default=0),
                "overlapping_fragment_window_observations": sum(int(r.get("overlapping_fragments", 0)) for r in sites),
                "discordant_fragment_window_observations": sum(int(r.get("discordant_fragments", 0)) for r in sites),
                "truth_snp_sites": len(truth_sites),
                "limited_truth_snp_sites": sum(r.get("cpg_context_limited", "false") == "true" for r in truth_sites),
                "truth_snp_sites_losing_all_links": sum(int(r["offered_linking_fragments"]) > 0 and int(r["retained_linking_fragments"]) == 0 for r in truth_sites),
                "missed_truth_alleles_with_limited_context": sum(r["label"] and float(r["posterior"]) < .5 and r.get("cpg_context_limited", "false") == "true" for r in data),
                "weighted_limited_fraction": sum(r["weight"] for r in limited) / site_weight if site_weight else None,
                "weighted_link_loss_fraction": lost_weight / offered_weight if offered_weight else None,
                "cpg_omitted_sites": sum(int(r.get("cpg_omitted_sites", 0)) for r in sites),
                "cpg_omitted_indels": sum(int(r.get("cpg_omitted_indels", 0)) for r in sites),
                "offered_links": offered, "retained_links": retained,
                "lost_links": offered - retained,
                "sites_losing_all_links": sum(int(r["offered_linking_fragments"]) > 0 and int(r["retained_linking_fragments"]) == 0 for r in sites),
                "context_changes_presence_calls": sum((float(r["posterior"]) >= .5) != (float(r["site_posterior"]) >= .5) for r in data),
                "alleles_with_probability_change_over_0_1": sum(abs(float(r["posterior"]) - float(r["site_posterior"])) > .1 for r in data),
                "max_probability_change": max((abs(float(r["posterior"]) - float(r["site_posterior"])) for r in data), default=0)}}


def evaluate_region(selected, manifest, root, binary, tabix, reuse=False):
    region = selected["region"]
    folder = root / region.split(":")[0]
    prepared = json.loads((folder / "prepared.json").read_text())
    if digest(folder / "sites.bed") != prepared["sites_sha256"]:
        raise ValueError("The predefined target panel changed")
    population = sum(r["weight"] for r in prepared["panel"])
    if abs(population - round(population)) > 1e-7:
        raise ValueError("Panel weights do not reconstruct an integer target population")
    prepared["scored_population_sites"] = round(population)
    inputs = manifest["inputs"]
    truth, gt, indels, _, _ = region_data(inputs, region, tabix)
    reference = reference_sequence(inputs["reference"], region)
    commands = json.loads((folder / "report.json").read_text())["commands"] if reuse else {}
    metrics = {}
    for name, library, engine in arms(manifest):
        output = folder / f"{name}.tsv"
        command = [str(binary), "ml", "haplotypes", "--fasta-file", inputs["reference"], inputs[library],
                   "--region", region, "--conversion-rate", ".95", "--false-conversion-rate", ".003",
                   "--min-alt-observations", "2", "--likelihood-engine", engine,
                   "--sites-bed", str(folder / "sites.bed"), "--output", str(output)]
        if library == "wgs":
            command.append("--non-taps")
        command.extend(arm_arguments(manifest, name))
        if not reuse:
            commands[name] = run(command, folder / name)
        data = rows(output)
        validate(data, {(r["chrom"], r["pos"]) for r in prepared["panel"]})
        metrics[name] = summarize(decorate(data, prepared, truth, gt, indels, reference))
        print(f"Finished {region} {name}: TP/FP/FN={tuple(metrics[name]['alleles'][k] for k in ['tp','fp','fn'])}", flush=True)
    gatk, _ = variants(query(tabix, inputs["gatk"], region))
    panel = {(r["chrom"], r["pos"]) for r in prepared["panel"]}
    gatk = {k for k in gatk if k[:2] in panel}
    report = {"prepared": prepared, "commands": commands, "metrics": metrics,
              "gatk_panel": {"tp": len(gatk & truth), "fp": len(gatk - truth), "fn": len(truth - gatk)}}
    (folder / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--rastair", type=pathlib.Path, required=True)
    parser.add_argument("--tabix", default="tabix")
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--score-only", action="store_true")
    parser.add_argument("--summarize-only", action="store_true", help="rescore existing TSVs without rerunning inference")
    parser.add_argument("--panel-source", type=pathlib.Path, help="reuse the exact target panels from an earlier suite")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    if args.jobs < 1 or manifest["negative_sites_per_region"] < 1 or sum([args.prepare_only, args.score_only, args.summarize_only]) > 1:
        parser.error("Invalid jobs, sample size or stage selection")
    chromosomes = [r["region"].split(":")[0] for r in manifest["regions"]]
    if len(set(chromosomes)) != len(chromosomes):
        parser.error("Use one interval per chromosome")
    root, binary = args.output_dir.resolve(), args.rastair.resolve()
    if args.panel_source and args.panel_source.resolve() == root:
        parser.error("Panel source and output directory must differ")
    root.mkdir(parents=True, exist_ok=True)
    metadata = {"manifest": manifest, "binary_sha256": digest(binary),
                "script_sha256": digest(__file__), "created_unix_time": time.time(),
                "probability_scope": "Candidate union plus truth SNPs and specified known disagreement sites in confident windows; negative candidate sites sampled with inverse-probability calibration weights.",
                "counts_scope": "TP/FP/FN and genotype errors refer to the evaluated panel. Estimated candidate-population FP counts are separately labelled.",
                "context_scope": "Lost links count fragment-site-window incidences, not unique molecules across the genome. Site-only ablation measures retained context, not an unrestricted multi-CpG oracle."}
    if not (args.score_only or args.summarize_only):
        (root / "suite.json").write_text(json.dumps(metadata, indent=2) + "\n")
        if args.panel_source:
            source = args.panel_source.resolve()
            original = json.loads((source / "suite.json").read_text())
            for field in ("inputs", "regions", "negative_sites_per_region", "seed", "known_sites"):
                if original["manifest"].get(field) != manifest.get(field):
                    raise ValueError(f"Panel selection differs in {field}")
            for chromosome in chromosomes:
                folder = root / chromosome
                folder.mkdir(exist_ok=True)
                prepared = json.loads((source / chromosome / "prepared.json").read_text())
                if digest(source / chromosome / "sites.bed") != prepared["sites_sha256"]:
                    raise ValueError("Source panel changed")
                prepared["panel_source"] = {"directory": str(source), "discovery_suite": original}
                shutil.copyfile(source / chromosome / "sites.bed", folder / "sites.bed")
                (folder / "prepared.json").write_text(json.dumps(prepared, indent=2) + "\n")
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
                futures = [pool.submit(prepare_region, s, manifest, root, binary, args.tabix) for s in manifest["regions"]]
                for future in futures:
                    future.result()
    else:
        old = json.loads((root / "suite.json").read_text())
        if old["manifest"] != manifest or old["binary_sha256"] != metadata["binary_sha256"]:
            raise ValueError("The binary or predefined manifest changed between discovery and scoring")
    if args.prepare_only:
        return
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        reports = list(pool.map(lambda s: evaluate_region(s, manifest, root, binary, args.tabix, args.summarize_only), manifest["regions"]))
    all_data = {name: [] for name, _, _ in arms(manifest)}
    for selected, report in zip(manifest["regions"], reports):
        truth, gt, indels, _, _ = region_data(manifest["inputs"], selected["region"], args.tabix)
        reference = reference_sequence(manifest["inputs"]["reference"], selected["region"])
        for name in all_data:
            all_data[name].extend(decorate(rows(root / selected["region"].split(":")[0] / f"{name}.tsv"), report["prepared"], truth, gt, indels, reference))
    combined = {"suite": metadata, "discovery_suite": json.loads((root / "suite.json").read_text()), "regions": reports,
                "combined": {name: summarize(data) for name, data in all_data.items()},
                "held_out": {name: summarize([r for r in data if r["split"] == "held_out"]) for name, data in all_data.items()}}
    (root / "report.json").write_text(json.dumps(combined, indent=2) + "\n")
    with (root / "summary.tsv").open("w") as handle:
        fields = ["arm", "tp", "fp", "fn", "genotype_errors", "genotyped_truth_sites", "weighted_log_loss", "weighted_brier", "limited_panel_sites", "weighted_limited_fraction", "weighted_link_loss_fraction"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for name, value in combined["combined"].items():
            writer.writerow({"arm": name, **{k: value["alleles"][k] for k in ["tp", "fp", "fn", "weighted_log_loss", "weighted_brier"]},
                "genotype_errors": value["genotypes"]["genotype_errors"], "genotyped_truth_sites": value["genotypes"]["truth_snp_sites"],
                "limited_panel_sites": value["context"]["limited_sites"],
                "weighted_limited_fraction": value["context"]["weighted_limited_fraction"],
                "weighted_link_loss_fraction": value["context"]["weighted_link_loss_fraction"]})
    with (root / "errors.tsv").open("w") as handle:
        fields = ["arm", "error", "chrom", "pos", "ref", "alt", "posterior", "site_posterior", "fragments", "observations", "context", "near_indel", "disconnected_sites", "resource_limited_sites", "pruned_haplotypes", "maximum_pruned_mass", "offered_linking_fragments", "retained_linking_fragments", "weight", "split"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for name, data in all_data.items():
            for row in data:
                called = float(row["posterior"]) >= .5
                if called != row["label"]:
                    writer.writerow({"arm": name, "error": "FP" if called else "FN", **{f: row.get(f, "") if f in ("disconnected_sites", "resource_limited_sites", "pruned_haplotypes", "maximum_pruned_mass") else row[f] for f in fields if f not in ("arm", "error")}})
    print(json.dumps({name: value["alleles"] for name, value in combined["combined"].items()}, indent=2))


if __name__ == "__main__":
    main()
