#!/usr/bin/env python3
"""Train matched chr19 forests and evaluate full-PairHMM filtering on held-out chr12."""
import argparse
import gzip
import hashlib
import json
import os
import platform
import subprocess
import time
from pathlib import Path


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def run(command, directory):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / 'command.json').write_text(json.dumps(command, indent=2) + '\n')
    started = time.time()
    with (directory / 'log.txt').open('w') as log:
        result = subprocess.run(
            ['/usr/bin/time', '-f', '%e\t%U\t%S\t%M', '-o', str(directory / 'resources.tsv'), *command],
            stdout=log, stderr=subprocess.STDOUT,
        )
    (directory / 'status.json').write_text(json.dumps(dict(
        completed=result.returncode == 0, returncode=result.returncode,
        started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(started)),
        wall_seconds=time.time() - started, slurm_job_id=os.environ.get('SLURM_JOB_ID'),
        hostname=platform.node(),
    ), indent=2) + '\n')
    result.check_returncode()


def require_complete(directory):
    if not json.loads((directory / 'status.json').read_text())['completed']:
        raise RuntimeError(f'Incomplete run: {directory}')


def intervals(path, chromosome, strict=False):
    opener = gzip.open if str(path).endswith('.gz') else open
    selected = []
    with opener(path, 'rt') as source:
        for line in source:
            fields = line.split()
            if not fields or fields[0].startswith(('#', 'track', 'browser')):
                continue
            if len(fields) < 3:
                raise ValueError(f'Malformed BED record: {line.rstrip()}')
            start, end = int(fields[1]), int(fields[2])
            if not 0 <= start < end:
                raise ValueError(f'Invalid BED coordinates: {line.rstrip()}')
            if fields[0] == chromosome:
                selected.append((start, end))
            elif strict:
                raise ValueError(f'{path} contains {fields[0]}; expected only {chromosome}')
    merged = []
    for start, end in sorted(selected):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def intersect_intervals(first, second):
    result = []
    a = b = 0
    while a < len(first) and b < len(second):
        start, end = max(first[a][0], second[b][0]), min(first[a][1], second[b][1])
        if start < end:
            result.append((start, end))
        if first[a][1] <= second[b][1]:
            a += 1
        else:
            b += 1
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['train', 'call', 'compare'])
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--rastair', type=Path, required=True)
    parser.add_argument('--baseline-model', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--threads', type=int, required=True)
    parser.add_argument('--aardvark', type=Path)
    parser.add_argument('--haplotype-error-profile', type=Path,
                        help='Calibrated JSON profile embedded in the haplotype forest')
    parser.add_argument('--training-targets', type=Path,
                        help='Optional chr19 BED restricting expensive feature extraction')
    parser.add_argument('--validation-targets', type=Path,
                        help='Optional chr12 BED restricting calling and accuracy assessment')
    parser.add_argument('--arms', nargs='+', choices=['bundled', 'standard', 'haplotypes'])
    parser.add_argument('--phase-calls', action='store_true',
                        help='Also emit final genotype-conditioned phase calls in every calling arm')
    args = parser.parse_args()
    for path, chromosome in [(args.training_targets, 'chr19'),
                             (args.validation_targets, 'chr12')]:
        if path is not None and not intervals(path, chromosome, strict=True):
            parser.error(f'{path} contains no {chromosome} targets')
    inputs = json.loads(args.inputs.read_text())
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.stage == 'compare':
        calling = json.loads((output / 'call-manifest.json').read_text())
        if calling.get('validation_targets'):
            target = Path(calling['validation_targets'])
            if digest(target) != calling['validation_targets_sha256']:
                raise ValueError('Validation targets changed after calling')
            if args.validation_targets is not None and digest(args.validation_targets) != digest(target):
                raise ValueError('Comparison targets differ from the called pilot targets')
            args.validation_targets = target
    binary = str(args.rastair.resolve())
    masks = [arg for path in inputs['phase_masks'] for arg in ['--phase-mask', path]]
    manifest = dict(inputs=inputs, binary=binary, binary_sha256=digest(binary),
                    baseline_model_sha256=digest(args.baseline_model), threads=args.threads,
                    stage=args.stage, slurm_job_id=os.environ.get('SLURM_JOB_ID'),
                    phase_calls=args.phase_calls,
                    engine='PairHMM at every eligible candidate; shared CpG shrinkage')
    if args.haplotype_error_profile is not None:
        manifest['error_profile'] = str(args.haplotype_error_profile.resolve())
        manifest['error_profile_sha256'] = digest(args.haplotype_error_profile)
    for name, path in [('training_targets', args.training_targets),
                       ('validation_targets', args.validation_targets)]:
        if path is not None:
            manifest[name] = str(path.resolve())
            manifest[name + '_sha256'] = digest(path)
    (output / f'{args.stage}-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    if args.stage == 'train':
        for arm in args.arms or ['standard', 'haplotypes']:
            if arm == 'bundled':
                parser.error('The bundled model is a comparison arm, not a training arm')
            directory = output / arm
            run([binary, 'ml', 'train', '-r', inputs['reference'], inputs['taps'],
                 inputs['training_truth'], '-R', inputs['training_regions'],
                 '-l', '@' + str(args.training_targets.resolve()) if args.training_targets else 'chr19',
                 '--ml-features', arm, '--seed', '20261014', '--threads', str(args.threads),
                 '--export-features', str(directory / 'features'),
                 '--feature-analytics', str(directory / 'analytics'),
                 '-o', str(directory / 'model.rff.mpk.lz4'),
                 *(['--haplotype-error-profile', str(args.haplotype_error_profile.resolve())]
                   if arm == 'haplotypes' and args.haplotype_error_profile is not None else []),
                 *masks], directory)
    elif args.stage == 'call':
        for arm, model in [('bundled', args.baseline_model.resolve()),
                           ('standard', output / 'standard/model.rff.mpk.lz4'),
                           ('haplotypes', output / 'haplotypes/model.rff.mpk.lz4')]:
            if args.arms is not None and arm not in args.arms:
                continue
            if arm != 'bundled':
                require_complete(output / arm)
            directory = output / 'chr12' / arm
            run([binary, 'call', '-r', inputs['reference'], inputs['taps'],
                 '-l', '@' + str(args.validation_targets.resolve()) if args.validation_targets else 'chr12',
                 '--model', str(model), '--experimental-indels', '--vcf-all-fields',
                 '--threads', str(args.threads), '--vcf', str(directory / 'calls.bcf'),
                 *(['--phase'] if args.phase_calls else []), *masks], directory)
    else:
        if args.aardvark is None:
            parser.error('--aardvark is required for compare')
        import pysam
        import pysam.bcftools
        bed = output / 'chr12/confident.bed'
        confident = intervals(inputs['training_regions'], 'chr12')
        if args.validation_targets is not None:
            confident = intersect_intervals(confident, intervals(args.validation_targets, 'chr12'))
        if not confident:
            raise ValueError('No confident chr12 evaluation intervals')
        bed.parent.mkdir(parents=True, exist_ok=True)
        bed.write_text(''.join(f'chr12\t{start}\t{end}\n' for start, end in confident))
        truth = output / 'chr12/truth.vcf.gz'
        pysam.bcftools.view('-r', 'chr12', inputs['training_truth'], '-Oz', '-o', str(truth), catch_stdout=False)
        pysam.tabix_index(str(truth), preset='vcf', csi=True, force=True)
        with pysam.VariantFile(str(truth)) as source:
            samples = list(source.header.samples)
        truth_sample = inputs.get('truth_sample')
        if truth_sample is None:
            if len(samples) != 1:
                raise ValueError('Set truth_sample in the inputs JSON for a multi-sample truth VCF')
            truth_sample = samples[0]
        if truth_sample not in samples:
            raise ValueError(f'Truth sample {truth_sample!r} is absent from {truth}')
        for arm in args.arms or ['bundled', 'standard', 'haplotypes']:
            directory = output / 'chr12' / arm
            require_complete(directory)
            query = directory / 'query.vcf.gz'
            selected = directory / 'pass.bcf'
            pysam.bcftools.index('-f', str(directory / 'calls.bcf'), catch_stdout=False)
            pysam.bcftools.view('-R', str(bed), '-f', 'PASS', '-i', 'GT="alt"', str(directory / 'calls.bcf'),
                               '-Ob', '-o', str(selected), catch_stdout=False)
            pysam.bcftools.norm('-f', inputs['reference'], '-m', '-any', str(selected),
                               '-Oz', '-o', str(query), catch_stdout=False)
            pysam.tabix_index(str(query), preset='vcf', csi=True, force=True)
            run([str(args.aardvark.resolve()), 'compare', '--reference', inputs['reference'],
                 '--truth-vcf', str(truth), '--truth-sample', truth_sample, '--query-vcf', str(query),
                 '--regions', str(bed), '--threads', str(args.threads),
                 '--output-dir', str(directory / 'accuracy')], directory / 'comparison')


if __name__ == '__main__':
    main()
