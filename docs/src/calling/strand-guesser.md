# Read Orientation Guessing (Strand Guesser)

By default, Rastair determines @OT/@OB orientation from SAM flags.
With `--guess-read-orientation`, Rastair instead infers OT/OB from mismatch motifs observed in each read.

## Why This Exists

Some datasets have read orientation metadata that is incomplete, inconsistent, or incompatible with expected TAPS pairing conventions.
This is especially common in tagmentation-based library preparations, which can produce non-directional reads.
The strand guesser provides an evidence-based fallback so methylation/variant strand assignments can still be made.

## Inference Algorithm

For each aligned read, Rastair scans mismatch positions against the reference.
At each mismatch, it inspects both 2 bp windows that include the mismatch base:

- current + next base
- previous + current base

Motifs are counted on the read sequence as represented by htslib/reference orientation:

- `TG` evidence supports OT
- `CA` evidence supports OB

Assignment rule:

- `TG > CA` => OT
- `CA > TG` => OB
- tie or no evidence => deterministic pseudo-random split per read

The tie-break is deterministic, so repeated runs produce stable OT/OB assignments for the same read.

## Scope

This option currently affects the main pileup-based `call` workflow.
