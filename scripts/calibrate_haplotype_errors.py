#!/usr/bin/env python3
"""Fit molecular repeat errors between intact flanks in homozygous-reference truth regions."""
import argparse
import bisect
import collections
import concurrent.futures
import gzip
import json
import math
import re
from pathlib import Path


def merge(intervals):
    result = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def overlaps(intervals, starts, start, end):
    i = bisect.bisect_left(starts, end) - 1
    return i >= 0 and intervals[i][1] > start


def converted_match(observed, expected, taps):
    # Either orientation may occur in a library. Requiring consistent conversion
    # across each flank avoids treating arbitrary substitutions as chemistry.
    if len(observed) != len(expected):
        return False
    return observed == expected or (taps and any(all(a == b or (a, b) == pair
        for a, b in zip(observed, expected)) for pair in [('T', 'C'), ('A', 'G')]))


def context_counts(task):
    import pysam
    bam_path, chromosome, reference, tracts, taps = task
    if not tracts:
        return {}
    starts = [tract[0] for tract in tracts]
    observations = collections.defaultdict(dict)
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam.fetch(chromosome, max(0, starts[0] - 10), tracts[-1][1] + 10):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail or not 30 <= read.mapping_quality < 255:
                continue
            if not read.query_name or read.query_sequence is None or read.query_qualities is None:
                continue
            lo = bisect.bisect_left(starts, read.reference_start + 10)
            hi = bisect.bisect_right(starts, read.reference_end - 10)
            offered = [(i, tract) for i, tract in enumerate(tracts[lo:hi], lo) if tract[1] + 10 <= read.reference_end]
            if not offered:
                continue
            wanted = {pos for _, (start, end, _) in offered for pos in (start - 10, end + 9)}
            positions = {r: q for q, r in read.get_aligned_pairs() if r in wanted and q is not None}
            identity = (read.get_tag('RG') if read.has_tag('RG') else '', read.query_name)
            for index, (start, end, period) in offered:
                left, right = positions.get(start - 10), positions.get(end + 9)
                if left is None or right is None or right - left < 19:
                    continue
                observed = read.query_sequence[left:right + 1]
                qualities = read.query_qualities[left:right + 1]
                if min(qualities) < 20 or max(qualities) == 255:
                    continue
                left_ref, right_ref = reference[start - 10:start], reference[end:end + 10]
                if not converted_match(observed[:10], left_ref, taps) or not converted_match(observed[-10:], right_ref, taps):
                    continue
                tract = observed[10:-10]
                motif = reference[start:start + period]
                expected = (motif * ((len(tract) + period - 1) // period))[:len(tract)]
                if not converted_match(tract, expected, taps):
                    continue
                delta = len(tract) - (end - start)
                # A duplicate fetch/secondary alignment cannot become another vote.
                molecule = observations[index].setdefault(identity, {})
                molecule[(read.reference_start, read.flag)] = delta
    totals = collections.defaultdict(lambda: [0.0] * len(COLUMNS))
    for index, molecules in observations.items():
        start, end, period = tracts[index]
        key = (period, min(20, (end - start) // period))
        for reads in molecules.values():
            values = list(reads.values())
            weight = 1.0 / len(values)
            totals[key][0] += end - start
            totals[key][1] += sum(weight for delta in values if delta > 0)
            totals[key][2] += sum(weight for delta in values if delta < 0)
            totals[key][3] += sum(weight * abs(delta) for delta in values if delta)
            totals[key][4] += 1
    return dict(totals)


COLUMNS = ['exposure', 'insertions', 'deletions', 'error_bases', 'molecules']


def fit_profile(counts, platform, library, prior_exposure=1000.0):
    repeats = []
    for period in range(1, 7):
        extension_events = sum(v[1] + v[2] for (p, _), v in counts.items() if p == period)
        extension_bases = sum(v[3] for (p, _), v in counts.items() if p == period)
        continuation = (max(0.0, extension_bases - extension_events) + 10.0) / (extension_bases + 100.0)
        extension = max(1, min(60, round(-10 * math.log10(continuation))))
        for copies in range(2, 21):
            baseline_q = max(10, round(40 - math.exp(copies / (3 * math.pi)) + 1))
            baseline = 10 ** (-baseline_q / 10)
            exposure, insertions, deletions = (list(counts.get((period, copies), [])) + [0.0] * 3)[:3]
            qualities = [max(6, min(60, round(-10 * math.log10((events + prior_exposure * baseline) / (exposure + prior_exposure)))))
                         for events in (insertions, deletions)]
            repeats.append(dict(period=period, minimum_copies=copies,
                insertion_open_phred=qualities[0], deletion_open_phred=qualities[1], extension_phred=extension))
    return dict(version=1, platform=platform, library=library, maximum_evidence_phred=45, repeats=repeats)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ['bam', 'reference', 'truth', 'confidence', 'output']:
        parser.add_argument('--' + arg, required=True)
    parser.add_argument('--chromosome', default='chr19')
    parser.add_argument('--platform', default='illumina')
    parser.add_argument('--non-taps', action='store_true')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--per-context', type=int, default=500)
    args = parser.parse_args()
    import pysam
    with pysam.FastaFile(args.reference) as fasta:
        reference = fasta.fetch(args.chromosome).upper()
    open_bed = gzip.open if args.confidence.endswith('.gz') else open
    with open_bed(args.confidence, 'rt') as handle:
        confidence = merge((int(fields[1]), int(fields[2])) for line in handle
            if (fields := line.split()) and fields[0] == args.chromosome)
    confidence_starts = [start for start, _ in confidence]
    with pysam.VariantFile(args.truth) as truth:
        excluded = merge((max(0, record.start - 20), record.stop + 20) for record in truth.fetch(args.chromosome))
    excluded_starts = [start for start, _ in excluded]
    candidates = []
    seen = set()
    strata = collections.Counter()
    for period in range(1, 7):
        pattern = re.compile(r'([ACGT]{' + str(period) + r'})\1{1,}')
        for match in pattern.finditer(reference):
            start, end = match.span()
            if end - start < 4 or (start, end) in seen:
                continue
            motif = match.group(1)
            if any(period % shorter == 0 and motif == motif[:shorter] * (period // shorter) for shorter in range(1, period)):
                continue
            i = bisect.bisect_right(confidence_starts, start - 10) - 1
            if i < 0 or confidence[i][1] < end + 10 or overlaps(excluded, excluded_starts, start - 10, end + 10):
                continue
            key = (period, min(20, (end - start) // period))
            # Evenly distributed deterministic selection rather than a prefix of chr19.
            digest = (start * 2654435761) & 0xffffffff
            candidates.append((digest, key, (start, end, period)))
            seen.add((start, end))
    selected = []
    for _, key, tract in sorted(candidates):
        if strata[key] < args.per_context:
            selected.append(tract)
            strata[key] += 1
    selected.sort()
    chunks = collections.defaultdict(list)
    for tract in selected:
        chunks[tract[0] // 1_000_000].append(tract)
    tasks = [(args.bam, args.chromosome, reference, tracts, not args.non_taps) for tracts in chunks.values()]
    totals = collections.defaultdict(lambda: [0.0] * len(COLUMNS))
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as pool:
        for finished, counts in enumerate(pool.map(context_counts, tasks), 1):
            for key, values in counts.items():
                totals[key] = [a + b for a, b in zip(totals[key], values)]
            print(f'Completed {finished}/{len(tasks)} calibration chunks', flush=True)
    profile = fit_profile(totals, args.platform, 'wgs' if args.non_taps else 'taps')
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(profile, indent=2) + '\n')
    output.with_suffix('.counts.json').write_text(json.dumps(dict(arguments=vars(args), tracts=len(selected),
        counts=[dict(period=k[0], copies=k[1], **dict(zip(COLUMNS, v)))
                for k, v in sorted(totals.items())]), indent=2) + '\n')


if __name__ == '__main__':
    main()
