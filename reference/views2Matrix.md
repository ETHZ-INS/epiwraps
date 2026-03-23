# views2Matrix

converts a RleViews or RleViewsList with views of the same width to a
matrix, setting out-of-bounds regions to NA (or \`padVal\`).

## Usage

``` r
views2Matrix(v, padVal = NA_integer_)
```

## Arguments

- v:

  A \`RleViews\` or \`RleViewsList\` object with views of the same
  width.

- padVal:

  The value to assign to out-of-bound regions.

## Value

A numeric matrix.

## Examples

``` r
# we create an example RleViews with out-of-bound regions:
library(IRanges)
co <- Rle(values=c(0,1,0), lengths=c(100,50,20))
v <- Views(co, c(25,150),c(50,175))
# convert to matrix
views2Matrix(v)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
#> [1,]    0    0    0    0    0    0    0    0    0     0     0     0     0     0
#> [2,]    1    0    0    0    0    0    0    0    0     0     0     0     0     0
#>      [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
#> [1,]     0     0     0     0     0     0     0     0     0     0     0     0
#> [2,]     0     0     0     0     0     0     0    NA    NA    NA    NA    NA
```
