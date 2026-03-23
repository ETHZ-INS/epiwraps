# showTrackInfo

Provide some information about the relative signal ranges of each track.

## Usage

``` r
showTrackInfo(x, assay = "input", doPrint = TRUE)
```

## Arguments

- x:

  A named list of signal matrices or an EnrichmentSE object as produced
  by
  [`signal2Matrix`](https://ethz-ins.github.io/epiwraps/reference/signal2Matrix.md)

- assay:

  The assay to use, defaults to the input assay.

- doPrint:

  Logical; whether to print the information.

## Value

An invisible list of captions.

## Examples

``` r
data(exampleESE)
showTrackInfo(exampleESE)
#> H3K27ac ( 150x80 ) :
#>   -2kb/+2kb (40 windows each)
#>   around the centers of given regions 
#> H3K4me3 ( 150x80 ) :
#>   -2kb/+2kb (40 windows each)
#>   around the centers of given regions 
#> p300 ( 150x80 ) :
#>   -2kb/+2kb (40 windows each)
#>   around the centers of given regions 
```
