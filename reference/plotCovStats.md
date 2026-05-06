# plotCovStats

Plots coverage statistics, such as as fingerprint plot.

## Usage

``` r
plotCovStats(qc, labels = "AUTO", show.legend = TRUE)
```

## Arguments

- qc:

  A list of coverage statistics, as produced by
  [`getCovStats`](https://ethz-ins.github.io/epiwraps/reference/getCovStats.md).

- labels:

  Passed to
  [`plot_grid`](https://wilkelab.org/cowplot/reference/plot_grid.html).

- show.legend:

  Logical; whether to show the plot legend.

## Value

A grid object to be plotted.

## Examples

``` r
# we use an example bigwig file
bwf <- system.file("extdata/example_atac.bw", package="epiwraps")
# because most of the file is empty, we'll exclude some of the ranges
cs <- getCovStats(bwf, exclude=GRanges("1", IRanges(1, 4300000)))
plotCovStats(cs)
```
