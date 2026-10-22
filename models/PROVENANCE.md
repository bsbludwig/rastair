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
  the second, cascade-shaped collection, whose negative pools are small (other 18.7k, cpg
  10.6k, denovo 3.0k), so their fits take seconds.
- Cost on an M4 Pro: 321 s wall, 2,248 CPU-s, 11.8 GiB peak RSS. File 30.3 MB.

The file was trained at `9d8392fa` (`feature/engine-beta`), whose SNV engine profiles the level
of a `CpG` a substitution creates or destroys and mixes each read with its mismapping chance, so
the `ev_*` features differ from the previous model's; `RASTAIR_ARTIFACT_FLOOR` unset. The engine
pools were other 18,673, cpg 10,559, denovo 3,035 negatives. md5 `876347408ca50443dfd073be7d9a52d4`.
