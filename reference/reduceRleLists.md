# reduceRleLists

reduceRleLists

## Usage

``` r
reduceRleLists(res, fn = "+")
```

## Arguments

- res:

  A list of RleList objects

- fn:

  The function to use to combine them

## Value

A combined RleList object.

## Examples

``` r
list_of_rlelists <- list(
  RleList(A=c(0,0,1,0,0), B=c(0,0,0,0,1)),
  RleList(A=c(1,1,1,0,0), C=c(0,1,1,1,1)) )
reduceRleLists(list_of_rlelists)
#> RleList of length 3
#> $A
#> numeric-Rle of length 5 with 3 runs
#>   Lengths: 2 1 2
#>   Values : 1 2 0
#> 
#> $B
#> numeric-Rle of length 5 with 2 runs
#>   Lengths: 4 1
#>   Values : 0 1
#> 
#> $C
#> numeric-Rle of length 5 with 2 runs
#>   Lengths: 1 4
#>   Values : 0 1
#> 
```
