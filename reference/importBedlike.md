# importBedlike

Imports a bed-like file as a GRanges object. Uses \`rtracklayer\` import
functions if possible, and falls back onto an import that's not
format-committed otherwise.

## Usage

``` r
importBedlike(x, ...)
```

## Arguments

- x:

  The path to a bed or bed-like file (can be gzipped)

- ...:

  passed to [`fread`](https://rdrr.io/pkg/data.table/man/fread.html)

## Value

A \`GRanges\` object

## Examples

``` r
# example bed file:
filepath <- system.file("extdata/example_peaks.bed", 
                        package="epiwraps")
b <- importBedlike(filepath)
```
