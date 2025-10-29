# Machine-learning models

Rastair supports using machine-learning models to improve the accuracy of methylation and variant calling.

This is a new feature and is still experimental.
It should work quite well on regular human DNA with 10-20x converage.

## Parameters

The machine-learning model gives a prediction for how likely a variant or methylation call is.
Using the `--ml` parameter, you can set a threshold for this prediction, e.g. `--ml 0.75`.

To use different model files, you can pass `--model-cpg <FILE>`, `--model-denovo-cpg <FILE>`, and `--model-other <FILE>` parameters.
Only one set of models is included with Rastair at the moment.

## Performance

Running machine-learning model is by far the slowest part of the Rastair pipeline.
Reducing the number of positions, e.g. only requesting @CpG sites (using `--cpg-only`), will help speed up the process.

## Training

For now, the machine-learning models are trained outside of Rastair.
In the future, Rastair might gain the ability to train models directly,
which makes it easier to adapt models to different use cases.

<!-- document expected accuracy -->
<!-- document using other models -->
