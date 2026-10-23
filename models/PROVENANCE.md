# Model provenance

`rastair_default.rff.mpk.lz4` is the model compiled into the binary (`include_bytes!` in
`src/call/ml.rs`): a `RastairFlatModel` with three forests (`cpg`, `denovo`, `other`), each with
its own Platt scaling, a feature-layout fingerprint, and a nested `snv_engine` set of the same
three forests that `call --snv-engine` re-scores the engine's SNV candidates with. Indels have no
model; the pair-HMM calls them. Inspect a file with
`cargo run --release --example model_stats -- <path>`.

## Recipe

A seqair release build, every flag at its default except those shown:

```
rastair ml train --snv-engine \
  -r hg38.fa.gz \
  aaf06a50d5be4e349ced8a29133b61a2.bam \
  HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  -R HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
  -l "chr1 chr6 chr11" --seed 20260820 -o models/rastair_default.rff.mpk.lz4 -@ 8
```

- **BAM**: `aaf06a50…`, a ~26x NA12878 TAPS library. A 73x-trained SNV model was worse at 26x
  (crossover near 45x), so the model is trained at the depth it is mostly used at. The second
  library, `cc655007…`, is a hold-out; chr12, chr19 and chr20 are out of sample.
- **Truth**: GIAB HG001 v4.2.1, with its confident BED as `-R`.
- **Defaults**: 800 trees, `--max-features 4`, `--max-depth 40`, `--min-samples-leaf 10`,
  8,000 positive / 200,000 negative examples per model.
- The top-level forests come from the engine-off collection. The engine forests are fitted on
  the second, cascade-shaped collection: every alt the first stage scores at
  `CASCADE_MIN_SCORE` or more, whatever its read count. Its negative pools are small (other
  25.6k, cpg 10.7k, denovo 7.5k), so their fits take seconds.
- Cost on an M4 Pro (`-@ 12`): 290 s wall, 2,663 CPU-s, 13.3 GiB peak RSS. File 32.0 MB.

The file was trained at `4c339a3b` (`feature/haplotypings`). Pools: top-level other
462,625 / 15,589,292, cpg 80,544 / 7,083,948, denovo 93,478 / 1,654,618 (positive /
negative); engine other 458,951 / 25,593, cpg 78,852 / 10,690, denovo 92,558 / 7,456. md5
`b99d5ade14744e4946a1d3d0c5172fde`.
