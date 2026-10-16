#!/usr/bin/env python3
"""Assess SNP phase against phased truth in confident, unmasked BED intervals."""
import argparse
import bisect
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else Path(path).open()


class Intervals:
    def __init__(self, paths):
        regions = defaultdict(list)
        for path in paths:
            with open_text(path) as source:
                for line in source:
                    fields = line.split()
                    if not fields or fields[0].startswith(('#', 'track', 'browser')):
                        continue
                    if fields[:3] == ['Chr', 'Start', 'End']:
                        continue
                    chrom, start, end = fields[:3]
                    start, end = int(start), int(end)
                    if not 0 <= start < end:
                        raise ValueError(f'Invalid BED interval: {line}')
                    regions[chrom].append((start, end))
        self.regions = {}
        for chrom, intervals in regions.items():
            merged = []
            for start, end in sorted(intervals):
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
                else:
                    merged.append((start, end))
            self.regions[chrom] = merged

    def contains(self, chrom, pos):
        regions = self.regions.get(chrom, [])
        index = bisect.bisect_right(regions, (pos, float('inf'))) - 1
        return index >= 0 and pos < regions[index][1]


def snps(path, chrom, end):
    with open_text(path) as source:
        for line in source:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[0] != chrom:
                continue
            pos = int(fields[1]) - 1
            if pos >= end:
                break
            if len(fields) < 10:
                raise ValueError('Incomplete VCF record')
            alleles = [fields[3], *fields[4].split(',')]
            values = dict(zip(fields[8].split(':'), fields[9].split(':')))
            gt = values.get('GT', '.')
            if '.' in gt:
                continue
            genotype = tuple(alleles[int(index)] for index in gt.replace('|', '/').split('/'))
            if len(genotype) == 2 and all(base in 'ACGT' and len(base) == 1 for base in genotype):
                yield pos, genotype, '|' in gt, values.get('PS', 'chromosome'), fields[6]


def phase_errors(blocks):
    pairs = switches = minority = assessable_sites = 0
    for block in blocks.values():
        orientations = [orientation for _, orientation in sorted(block)]
        if len(orientations) < 2:
            continue
        assessable_sites += len(orientations)
        pairs += len(orientations) - 1
        switches += sum(a != b for a, b in zip(orientations, orientations[1:]))
        counts = Counter(orientations)
        minority += min(counts[0], counts[1])
    return dict(assessable_sites=assessable_sites, adjacent_pairs=pairs,
                switch_errors=switches, switch_error_rate=switches / pairs if pairs else None,
                minimum_orientation_disagreements=minority,
                orientation_disagreement_rate=minority / assessable_sites if assessable_sites else None,
                assessable_blocks=sum(len(block) >= 2 for block in blocks.values()))


def assess(truth_path, query_path, chrom, end, confidence, masks):
    eligible = lambda pos: confidence.contains(chrom, pos) and not masks.contains(chrom, pos)
    truth = {pos: (gt, ps) for pos, gt, phased, ps, filt in snps(truth_path, chrom, end)
             if eligible(pos) and phased and gt[0] != gt[1] and filt in ('PASS', '.')}
    counts = Counter()
    matched = set()
    blocks = defaultdict(list)
    for pos, gt, phased, ps, filt in snps(query_path, chrom, end):
        if not eligible(pos):
            continue
        if gt[0] != gt[1]:
            counts['query_heterozygous_snps'] += 1
        if phased:
            counts['query_phased_snps'] += 1
        if pos not in truth:
            if phased:
                counts['phased_sites_without_heterozygous_truth'] += 1
            continue
        expected, truth_ps = truth[pos]
        if sorted(gt) != sorted(expected):
            counts['genotype_disagreements_at_truth_hets'] += 1
            if phased:
                counts['phased_genotype_disagreements'] += 1
            continue
        matched.add(pos)
        if phased:
            counts['matching_phased_truth_hets'] += 1
            blocks[(ps, truth_ps)].append((pos, int(gt[0] != expected[0])))
    return dict(truth=str(truth_path), query=str(query_path), region=f'{chrom}:1-{end}',
                confident_unmasked_truth_hets=len(truth), matching_called_truth_hets=len(matched),
                phase_coverage_of_truth=counts['matching_phased_truth_hets'] / len(truth) if truth else None,
                phase_coverage_of_matching_calls=counts['matching_phased_truth_hets'] / len(matched) if matched else None,
                **counts, **phase_errors(blocks))


def compare_shared(truth_path, query_paths, chrom, end, confidence, masks):
    truth = {pos: (gt, ps) for pos, gt, phased, ps, filt in snps(truth_path, chrom, end)
             if phased and gt[0] != gt[1] and filt in ('PASS', '.')
             and confidence.contains(chrom, pos) and not masks.contains(chrom, pos)}
    queries = [{pos: (gt, ps) for pos, gt, phased, ps, filt in snps(path, chrom, end)
                if phased and pos in truth and sorted(gt) == sorted(truth[pos][0])}
               for path in query_paths]
    shared = set.intersection(*(set(query) for query in queries))
    arms = {}
    for path, query in zip(query_paths, queries):
        blocks = defaultdict(list)
        for pos in shared:
            key = (*[other[pos][1] for other in queries], truth[pos][1])
            blocks[key].append((pos, int(query[pos][0][0] != truth[pos][0][0])))
        arms[str(path)] = phase_errors(blocks)
    return dict(shared_genotype_matching_phased_sites=len(shared),
                definition='All arms assessed on the same adjacent SNP pairs within the intersection of their phase sets and truth phase sets.',
                arms=arms)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--truth', type=Path, required=True)
    parser.add_argument('--query', type=Path, nargs='+', required=True)
    parser.add_argument('--confidence-bed', type=Path, required=True)
    parser.add_argument('--exclude-bed', type=Path, action='append', default=[])
    parser.add_argument('--chromosome', default='chr1')
    parser.add_argument('--end', type=int, required=True, help='Exclusive zero-based end (equal to last included one-based base)')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    confidence, masks = Intervals([args.confidence_bed]), Intervals(args.exclude_bed)
    report = dict(confidence_bed=str(args.confidence_bed), exclude_bed=list(map(str, args.exclude_bed)),
                  definition='Switch errors are orientation changes between adjacent genotype-matched phased SNPs in the same query and truth phase set; isolated flips count as two switches. Genotype errors are excluded and reported separately.',
                  runs=[assess(args.truth, path, args.chromosome, args.end, confidence, masks) for path in args.query])
    if len(args.query) > 1:
        report['matched_link_comparison'] = compare_shared(args.truth, args.query, args.chromosome, args.end, confidence, masks)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
