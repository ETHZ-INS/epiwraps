# resize a numeric matrix to given dimensions

resize a numeric matrix to given dimensions

## Usage

``` r
resizeMatrix(mat, ndim = dim(mat), method = c("mean", "max", "min"))
```

## Arguments

- mat:

  A numeric matrix

- ndim:

  The desired output dimensions

- method:

  Whether to use normal interpolation (\`method="mean"\`, the default),
  or the max or min of the overlapping grid cells.

## Value

A numeric matrix of dimensions \`ndim\`

## Details

For most cased this is based on Vyha's implementation (taken from
https://stackoverflow.com/a/23429527 ), but adding the possibility to
replace normal interpolation with the max/min of the overlapping grid
cells. This works well when the desired dimensions are at least half of
the input ones. When the desired dimensions are smaller, a different
binning method is used, first applied on columns and then on rows.
