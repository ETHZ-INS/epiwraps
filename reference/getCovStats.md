# getCovStats

Assembles read distribution statistics from a set of bigwig files based
on random windows.

## Usage

``` r
getCovStats(
  x,
  binSize = 1000,
  nbBins = 10000,
  exclude = NULL,
  canonical.chr = TRUE,
  maxCovQuant = 0.999,
  BPPARAM = SerialParam()
)
```

## Arguments

- x:

  A (named) vector of paths to bigwig files (all from the same genome)

- binSize:

  The size of bins

- nbBins:

  The approximate number of random bins. More bins gives more accurate
  readouts but take longer to read and compute.

- exclude:

  Region to exclude

- canonical.chr:

  Logical; whether to restrict the sampling to standard chromosomes.

- maxCovQuant:

  The quantile to use as maximum coverage (default 0.999)

- BPPARAM:

  BioParallel BPPARAM for multithreading across files.

## Value

A named list of QC tables

## Examples

``` r
# we use an example bigwig file
bwf <- system.file("extdata/example_atac.bw", package="epiwraps")
# because most of the file is empty, we'll exclude some of the ranges
cs <- getCovStats(bwf, exclude=GRanges("1", IRanges(1, 4300000)))
plotCovStats(cs)
```
