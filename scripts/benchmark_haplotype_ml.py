#!/usr/bin/env python3
"""Train matched chr19 forests and evaluate full-PairHMM filtering on held-out chr12."""
import argparse
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['train', 'call', 'compare'])
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--rastair', type=Path, required=True)
    parser.add_argument('--baseline-model', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--threads', type=int, required=True)
    parser.add_argument('--aardvark', type=Path)
    parser.add_argument('--arms', nargs='+', choices=['bundled', 'canonical', 'haplotypes'])
    args = parser.parse_args()
    inputs = json.loads(args.inputs.read_text())
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    binary = str(args.rastair.resolve())
    masks = [arg for path in inputs['phase_masks'] for arg in ['--phase-mask', path]]
    manifest = dict(inputs=inputs, binary=binary, binary_sha256=digest(binary),
                    baseline_model_sha256=digest(args.baseline_model), threads=args.threads,
                    stage=args.stage, slurm_job_id=os.environ.get('SLURM_JOB_ID'),
                    engine='PairHMM at every eligible candidate; shared CpG shrinkage')
    (output / f'{args.stage}-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    if args.stage == 'train':
        for arm in args.arms or ['canonical', 'haplotypes']:
            if arm == 'bundled':
                parser.error('The bundled model is a comparison arm, not a training arm')
            directory = output / arm
            run([binary, 'ml', 'train', '-r', inputs['reference'], inputs['taps'],
                 inputs['training_truth'], '-R', inputs['training_regions'], '-l', 'chr19',
                 '--ml-features', arm, '--seed', '20261014', '--threads', str(args.threads),
                 '--export-features', str(directory / 'features'),
                 '--feature-analytics', str(directory / 'analytics'),
                 '-o', str(directory / 'model.rff.mpk.lz4'), *masks], directory)
    elif args.stage == 'call':
        for arm, model in [('bundled', args.baseline_model.resolve()),
                           ('canonical', output / 'canonical/model.rff.mpk.lz4'),
                           ('haplotypes', output / 'haplotypes/model.rff.mpk.lz4')]:
            if args.arms is not None and arm not in args.arms:
                continue
            if arm != 'bundled':
                require_complete(output / arm)
            directory = output / 'chr12' / arm
            run([binary, 'call', '-r', inputs['reference'], inputs['taps'], '-l', 'chr12',
                 '--model', str(model), '--experimental-indels=ml', '--vcf-all-fields',
                 '--threads', str(args.threads), '--vcf', str(directory / 'calls.bcf'),
                 *masks], directory)
    else:
        if args.aardvark is None:
            parser.error('--aardvark is required for compare')
        import pysam
        import pysam.bcftools
        bed = output / 'chr12/confident.bed'
        bed.write_text(''.join(line for line in Path(inputs['training_regions']).read_text().splitlines(True)
                               if line.split('\t', 1)[0] == 'chr12'))
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
        for arm in args.arms or ['bundled', 'canonical', 'haplotypes']:
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
