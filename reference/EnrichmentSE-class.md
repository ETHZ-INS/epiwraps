# EnrichmentSE class and constructor

The \`EnrichmentSE\` class is a container for epigenomic enrichment
data, extending the
[`RangedSummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
class.

## Usage

``` r
EnrichmentSE(assays, rowRanges = NULL, ...)

# S4 method for class 'EnrichmentSE,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'EnrichmentSE'
show(object)

# S4 method for class 'EnrichmentSE'
score(x, ...)
```

## Arguments

- assays:

  A list of matrices or \`normalizedMatrix\` objects.

- rowRanges:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object.

- ...:

  Arguments passed to
  [`getSignalMatrices`](https://ethz-ins.github.io/epiwraps/reference/getSignalMatrices.md).

- x:

  An object of class \`EnrichmentSE\`, as produced by
  [`signal2Matrix`](https://ethz-ins.github.io/epiwraps/reference/signal2Matrix.md).

- i, j:

  Indices for subsetting.

- drop:

  Logical; whether to drop dimensions.

- object:

  An \`EnrichmentSE\` object.

## Value

An \`EnrichmentSE\` object.

## Functions

- `x[i`: Subset method for EnrichmentSE

- `show(EnrichmentSE)`: Show method for EnrichmentSE

- `score(EnrichmentSE)`: Access the score (first assay)
