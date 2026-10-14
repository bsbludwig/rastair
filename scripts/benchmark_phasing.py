#!/usr/bin/env python3
"""Time matched call/phase runs and verify that phasing preserves all other VCF fields."""
import argparse
import hashlib
import json
import platform
import subprocess
from collections import Counter
from pathlib import Path


def inspect_vcf(path):
    digest = hashlib.sha256()
    blocks = Counter()
    anchors = {}
    hets = phased = records = 0
    last = None
    with path.open() as source:
        for line in source:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            records += 1
            last = [fields[0], int(fields[1])]
            names = fields[8].split(':')
            values = dict(zip(names, fields[9].split(':')))
            gt = values.get('GT', '.')
            alleles = gt.replace('|', '/').split('/')
            het = len(alleles) == 2 and '.' not in alleles and alleles[0] != alleles[1]
            hets += het
            if '|' in gt:
                sequences = [fields[3], *fields[4].split(',')]
                if not het or any(len(sequences[int(allele)]) != 1 for allele in alleles):
                    raise ValueError(f'Phased non-SNP heterozygote: {line}')
                ps = int(values['PS'])
                key = (fields[0], ps)
                anchors.setdefault(key, (int(fields[1]), alleles))
                blocks[key] += 1
                phased += 1
            values['GT'] = '/'.join(sorted(alleles))
            values.pop('PS', None)
            digest.update(json.dumps([fields[:8], sorted(values.items())], separators=(',', ':')).encode())
    for key, size in blocks.items():
        pos, alleles = anchors[key]
        if size < 2 or pos != key[1] or int(alleles[0]) >= int(alleles[1]):
            raise ValueError(f'Invalid phase block: {key}')
    return dict(records=records, heterozygotes=hets, phased_sites=phased, phase_blocks=len(blocks), largest_block=max(blocks.values(), default=0), last_record=last, nonphase_sha256=digest.hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--rastair', type=Path, required=True)
    parser.add_argument('--legacy', type=Path)
    parser.add_argument('--regions', nargs='+', required=True)
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--phase-mask', type=Path, action='append', default=[])
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    inputs = json.loads(args.inputs.read_text())
    masks = args.phase_mask or list(map(Path, inputs.get('phase_masks', [])))
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    cpu_model = next((line.split(':', 1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')), platform.processor())
    report = dict(inputs=inputs, threads=args.threads, platform=platform.platform(), cpu_model=cpu_model, cache_policy='normal OS cache; unphased arm precedes phase arms', timing_scope='caller only, including input and VCF output; excludes subsequent VCF inspection', runs=[])
    arms = [('unphased', args.rastair, False), ('pair_hmm', args.rastair, True)]
    if args.legacy:
        arms.append(('legacy', args.legacy, True))
    for region in args.regions:
        region_runs = []
        for arm, binary, phase in arms:
            stem = region.replace(':', '_').replace('-', '_') + '_' + arm
            output = out / (stem + '.vcf')
            resource = out / (stem + '.resource.tsv')
            command = [str(binary.resolve()), 'call', inputs['taps'], '-r', inputs['reference'], '-l', region, '--threads', str(args.threads), '--vcf', str(output)]
            if phase:
                command.append('--phase')
                if arm != 'legacy':
                    for mask in masks:
                        command.extend(['--phase-mask', str(mask.resolve())])
            timed = ['/usr/bin/time', '-f', '%e\t%U\t%S\t%M', '-o', str(resource), *command]
            print(f'Running {arm}: {region}', flush=True)
            with (out / (stem + '.log')).open('w') as log:
                subprocess.run(timed, stdout=log, stderr=log, check=True)
            wall, user, system, rss = resource.read_text().split()
            run = dict(region=region, arm=arm, command=command, binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(), wall_seconds=float(wall), cpu_seconds=float(user)+float(system), peak_rss_kib=int(rss), **inspect_vcf(output))
            region_runs.append(run)
            report['runs'].append(run)
            (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
            print(json.dumps(run), flush=True)
        if len({run['nonphase_sha256'] for run in region_runs}) != 1:
            raise ValueError(f'Phasing changed nonphase fields in {region}')
    report['nonphase_fields_identical'] = True
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
