ResMDS
======

MDS uniform re-sampler for continuously differentiable block distance matrices.

use as a function [Y,stress]=resmds(M, d, resample, structure, showresult, embeddingfun)

[Y,stress]=resmds(M, d, resample, structure)

M - similarity matrix.

d - target space dimension.

resample - resampling ratio (5 -> 1/5 of M's size)

structure - Array of block sizes in M (each value corresponds to the dimensions of a block sub-matrix.

showresult - default is false. uses a scatter plot to show the first three-dimensions of resulting embedding.

embeddingfun - default is @mdscale. Holds a function handle to the embedding function.

