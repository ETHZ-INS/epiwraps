# computeDeviationsFaster

A faster version of `computeDeviations` (see details). This is kept for
legacy/backward compatibility, and users should rather use the
\`betterChromVAR\` package.

## Usage

``` r
computeDeviationsFaster(
  counts,
  motifMatches,
  backgrounds,
  normalize = TRUE,
  welford = NULL,
  verbose = TRUE,
  BPPARAM = SerialParam(progress = TRUE)
)
```

## Arguments

- counts:

  A matrix of read counts per region(rows)/sample(columns), or a
  SummarizedExperiment with this as first assay, as produced by
  `getCounts` or
  [`peakPbCountsSE`](https://ethz-ins.github.io/epiwraps/reference/peakPbCountsSE.md).
  Can also be already normalized if \`normalize=FALSE\`.

- motifMatches:

  A matrix of motif matches (either logical or with something akin to
  binding probabilities), or a SummarizedExperiment containing such an
  assay.

- backgrounds:

  A matrix of indices indicating background peaks, as produced by
  `getBackgroundPeaks`.

- normalize:

  Logical; whether to perform column sum normalization on \`counts\`

- welford:

  Logical; whether to use Welford's online algorithm for computing
  background means and standard deviations. If NULL (default), the
  function will use Welford's if the predicted memory usage is above
  10GB.

- verbose:

  Logical; whether to print messages, in particular the projected memory
  size for large datasets.

- BPPARAM:

  An optional BiocParallel BPPARAM object for multi-threading.

## Value

A SummarizedExperiment

## Details

The results of this function are equivalent (within a very small
precision) to those of `computeDeviations`, but faster due to its
reliance on matrix operations, which enables the use of a larger number
of background iterations. The use of matrix operations however implies
that the entire results of each background iteration is stored in
memory. While this is typically not a problem for bulk data, it often is
with large (e.g. single-cell) datasets. For this reason, the function
includes a variant using Welford's algorithm for the computation of the
means and standard deviations. This is still faster than the original
chromVAR implementation, but much less so than when using
\`welford=FALSE\`. By default, Welford's algorithm will be used if the
projected memory usage is above 10GB. The speed gain depends a lot on
the dimensions of the different inputs and the number of threads used.
The runtime when using Welford's algorithm can vary from 50 to 100
Welford's algorithm, running times will typically be much smaller
(typically roughly 20

## Author

Pierre-Luc Germain
