# regionOverlaps

A wrapper for visualizing pairwise-wise overlaps across multiple sets of
genomic ranges.

## Usage

``` r
regionOverlaps(
  listOfRegions,
  ignore.strand = TRUE,
  cluster = length(listOfRegions) > 2,
  color = viridis::plasma(100),
  number_color = "black",
  ...
)
```

## Arguments

- listOfRegions:

  A named list of two or more (non-empty) \`GRanges\`

- ignore.strand:

  Logical; whether to ignore strand for overlaps

- cluster:

  Logical; whether to cluster rows/columns

- color:

  Heatmap colorscale

- number_color:

  Values color

- ...:

  Passed to
  [`pheatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/pheatmap.html)

## Value

A \`Heatmap\` showing the overlap coefficient as colors, and the overlap
size as values.

## Examples

``` r
# random list of GRanges:
grl <- lapply(c(A=10,B=20,C=30), FUN=function(x){
  GRanges("seq1", IRanges(runif(x,1,1000), width=20))
})
regionOverlaps(grl)
#> Error in loadNamespace(x): there is no package called 'viridis'
```
