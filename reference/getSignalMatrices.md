# getSignalMatrices

Extracts a list of signal matrices from an EnrichmentSE object.

## Usage

``` r
getSignalMatrices(x, assay = 1L)
```

## Arguments

- x:

  An object of class \`EnrichmentSE\`, as produced by
  [`signal2Matrix`](https://ethz-ins.github.io/epiwraps/reference/signal2Matrix.md).

- assay:

  The assay to extract (defaults to the first assay).

## Value

A list of normalizedMatrix objects.

## Examples

``` r
# we first get an EnrichmentSE object:
data(exampleESE)
# then we can extract the list of signal matrices:
sm <- getSignalMatrices(exampleESE)
```
