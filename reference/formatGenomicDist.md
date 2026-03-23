# formatGenomicDist

formatGenomicDist

## Usage

``` r
formatGenomicDist(
  e,
  allowFraction = TRUE,
  sameUnits = TRUE,
  head0 = TRUE,
  units = c(mb = 1e+06, kb = 1000, bp = 1)
)
```

## Arguments

- e:

  An integer vector of genomic sizes/distances to be formatted

- allowFraction:

  Whether to allow decimals; can be either logical or a number of
  decimals allowed (defaults to TRUE / 1)

- sameUnits:

  Logical; whether all values should get the same unit.

- head0:

  Logical; whether to keep heading zero

- units:

  The units to use.

## Value

A character vector.
