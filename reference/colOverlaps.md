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
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    1
#> [3,]    0    1    1
```
