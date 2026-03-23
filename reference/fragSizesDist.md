# fragSizesDist

fragSizesDist

## Usage

``` r
fragSizesDist(
  x,
  what = 10000,
  flags = scanBamFlag(isProperPair = TRUE),
  BPPARAM = SerialParam(),
  returnPlot = TRUE
)
```

## Arguments

- x:

  A (named) vector of paths to bam files.

- what:

  Either a positive integer (length 1) indicating how many reads to
  randomly sample, or a character vector (of length 1) indicating which
  chromosome to read.

- flags:

  A \`scanBamFlag\` object (see
  [ScanBamParam](https://rdrr.io/pkg/Rsamtools/man/ScanBamParam-class.html))

- BPPARAM:

  A
  [BiocParallel](https://rdrr.io/pkg/BiocParallel/man/BiocParallel-package.html)
  BPPARAM object for multithreading.

- returnPlot:

  Logical; whether to return a plot.

## Value

If \`returnPlot=TRUE\`, returns a ggplot object, otherwise a data.frame
of fragment lengths.

## Examples

``` r
# example bam file:
bam <- system.file("extdata", "ex1.bam", package="Rsamtools")
fragSizesDist(bam, what=100)
```
