# Model provenance

`rastair_default.rff.mpk.lz4` is the model compiled into the binary (`include_bytes!` in
`src/call/ml.rs`): a `RastairFlatModel` holding one feature-layout fingerprint and two
`ForestSet`s of three forests (`cpg`, `denovo`, `other`), each with its own Platt scaling:
`forests`, which score every alt, and `snv_engine`, which `call --snv-engine` re-scores the
engine's SNV candidates with. Indels have no model; the pair-HMM calls them.

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

- **BAM**: `aaf06a50…`, a ~26x NA12878 TAPS library; its chr1, chr6 and chr11 are also in
  `bench/harness/data/aaf_train_26x.bam` (rastair2-merge-bench), which trains the same model. A 73x-trained SNV model was worse at 26x
  (crossover near 45x), so the model is trained at the depth it is mostly used at. The second
  library, `cc655007…`, is a hold-out; chr12, chr19 and chr20 are out of sample.
- **Truth**: GIAB HG001 v4.2.1, with its confident BED as `-R`; `ml train` needs the VCF's
  `.csi` (`tmp/1000genomes/` has one).
- **Defaults**: 800 trees, `--max-features 4`, `--max-depth 40`, `--min-samples-leaf 10`,
  8,000 positive / 200,000 negative examples per model.
- The top-level forests come from the engine-off collection. The engine forests are fitted on
  the second, cascade-shaped collection: every alt the first stage scores at
  `CASCADE_MIN_SCORE` or more, whatever its read count. Its negative pools are small (other
  23.4k, cpg 10.9k, denovo 6.6k), so their fits take seconds.
- Both collections take a candidate's row through `call::process::ml_row`, the gate `call`
  scores through.
- Cost on an M4 Pro (`-@ 8`): 406 s wall, 2,341 CPU-s user, 8.3 GiB peak RSS. File 31.2 MB.

The file was trained at `1090d5fa` from `bench/harness/data/aaf_train_26x.bam`, after the
H2/M8 read counts behind `ev_ad_*` and the M5 neighbour features changed. Pools: top-level
other 462,242 / 15,586,896, cpg 80,453 / 7,080,126, denovo 93,436 / 1,654,438 (positive /
negative); engine other 458,767 / 23,409, cpg 78,897 / 10,911, denovo 92,540 / 6,573. md5
`2784ef8b34fa726a118b0120eaf18372`.
