# Compute resources

## CPU

Run time drops almost linearly with `--threads` up to the number of physical cores.
We recommend using at most one thread per hyperthreading core.

## Memory

Rastair's peak memory grows with the number of threads,
the segment length (`--segment-max-length`, default 100 kb),
and the coverage:

| threads | 30×, 100 kb | 60×, 100 kb | 30×, 400 kb |
| ------: | ----------: | ----------: | ----------: |
|       8 |      0.7 GB |      0.9 GB |      1.3 GB |
|      16 |      1.1 GB |      1.4 GB |      2.3 GB |
|      32 |      1.9 GB |      2.6 GB |      4.2 GB |

```admonish note
These numbers are for the experimental seqair backend.
The default (htslib) backend needs about twice as much per thread, independent of coverage.
```

To reduce memory usage, lower `--threads` or `--segment-max-length`.

## GPU

`--gpu` speeds up the @ML predictions.
A whole run gets about 2x faster on an M4 MacBook Pro and about 2.2× faster on Linux with an AMD Radeon RX 5700 XT (16 threads).
It needs little extra RAM: none on macOS, about 130 MB on Linux.
See [GPU acceleration](./gpu.md) for the requirements.
