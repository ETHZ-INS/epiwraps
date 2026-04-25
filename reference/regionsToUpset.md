# regionsToUpset

Prepares sets of regions for UpSet overlap representation. A wrapper
around `upset` for comparing multiple sets of genomic ranges.

## Usage

``` r
regionsToUpset(
  x,
  reference = c("reduce", "disjoin"),
  returnList = FALSE,
  ignore.strand = FALSE,
  maxgap = -1L,
  minoverlap = 0L,
  ...
)
```

## Arguments

- x:

  A named list of genomic ranges (or paths to bed files)

- reference:

  The method for creating the reference windows ('reduce' or 'disjoin').
  Alternatively, a \`GRanges\` object of reference windows.

- returnList:

  Logical; whether to return the list of regions instead of plotting.

- ...:

  Further arguments specifying how the overlaps are done, passed to
  [`findOverlaps-methods`](https://rdrr.io/pkg/GenomicRanges/man/findOverlaps-methods.html)).

## Value

A data.frame of set inclusions which can be directly input to
[`make_comb_mat`](https://rdrr.io/pkg/ComplexHeatmap/man/make_comb_mat.html),
and then [`UpSet`](https://rdrr.io/pkg/ComplexHeatmap/man/UpSet.html).

## Examples

``` r
# random list of GRanges:
grl <- lapply(c(A=10,B=20,C=30), FUN=function(x){
  GRanges("seq1", IRanges(runif(x,1,1000), width=20))
})
input_for_upset <- regionsToUpset(grl)
# we would then plot the data with:
ComplexHeatmap::UpSet(make_comb_mat(input_for_upset))
```
