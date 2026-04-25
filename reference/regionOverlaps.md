# regionOverlaps

A wrapper for visualizing pairwise-wise overlaps across multiple sets of
genomic ranges.

## Usage

``` r
regionOverlaps(
  listOfRegions,
  mode = c("reduced", "pairwise"),
  ignore.strand = TRUE,
  cluster = length(listOfRegions) > 2,
  colorBy = c("overlapCoef", "jaccard"),
  ...,
  returnValues = FALSE,
  color = viridisLite::plasma(100),
  number_color = "black"
)
```

## Arguments

- listOfRegions:

  A named list of two or more (non-empty) \`GRanges\`

- mode:

  Either 'reduced' or 'pairwise'. 'reduced' first uses \`reduce\` to get
  a set of reference regions which are, based on overlap, contained or
  not in the different sets. It is thus symmetrical. \`pairwise\` does
  pairwise overlap between the sets of regions; it is asymmetrical and
  slower to compute.

- ignore.strand:

  Logical; whether to ignore strand for overlaps (default TRUE).

- cluster:

  Logical; whether to cluster rows/columns

- colorBy:

  Whether to color by 'overlapCoef' (default), or by 'jaccard' index.

- ...:

  Passed to
  [`pheatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/pheatmap.html)

- returnValues:

  Logical; whether to return the matrix of requested values instead of
  plotting.

- color:

  Heatmap colorscale

- number_color:

  Values color

## Value

A \`Heatmap\` showing the overlap coefficient as colors, and the overlap
size as values. Alternatively, if \`returnValues=TRUE\`, a matrix of the
requested values.

## Examples

``` r
# random list of GRanges:
grl <- lapply(c(A=10,B=20,C=30), FUN=function(x){
  GRanges("seq1", IRanges(runif(x,1,1000), width=20))
})
regionOverlaps(grl)
```
