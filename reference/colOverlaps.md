# colOverlaps

Computes pairwise overlap metric between columns

## Usage

``` r
colOverlaps(x, y = NULL, metric = c("overlap", "jaccard", "overlapCoef"))
```

## Arguments

- x:

  A logical matrix.

- y:

  An optional nother logical matrix or vector. If NULL (default),
  overlaps are computed between columns of \`x\`.

- metric:

  Either 'overlap', 'jaccard', or 'overlapCoef'.

## Value

A matrix of pairwise metric values.

## Examples

``` r
m <- matrix(sample(c(TRUE,FALSE),12,replace=TRUE), nrow=4)
colOverlaps(m, metric="jaccard")
#> Error in colSums(x): 'x' must be an array of at least two dimensions
```
