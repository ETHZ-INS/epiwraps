# meltSignals

Aggregates and melts a list of signal matrices, for plotting (with
ggplot).

## Usage

``` r
meltSignals(ml, fun = NULL, splitBy = NULL, trim = 0.98, assay = 1L)
```

## Arguments

- ml:

  A named list of signal matrices or an EnrichmentSE object as produced
  by
  [`signal2Matrix`](https://ethz-ins.github.io/epiwraps/reference/signal2Matrix.md)

- fun:

  An optional custom aggregation function (or named list thereof).

- splitBy:

  A vector of values (factor or character of length equal to
  \`nrow(ml)\`) by which to split the aggregation. Can also be the name
  of a column of \`rowData(ml)\`.

- trim:

  The quantile above which to trim values. If a numeric vector of length
  2, will be used as lower and upper quantiles beyond which to trim.

- assay:

  Assay to use (ignored unless \`ml\` is an ESE object), defaults to the
  first assay.

## Value

A data.frame.

## Examples

``` r
# we first get an EnrichmentSE object:
data(exampleESE)
# we extract the means per position:
d <- meltSignals(exampleESE)
head(d)
#>   position  sample     mean       SD       SE median
#> 1    -2000 H3K27ac 83.28000 182.5112 14.90197   11.5
#> 2    -1950 H3K27ac 85.87333 189.9628 15.51040   11.5
#> 3    -1900 H3K27ac 91.19333 207.2911 16.92525    9.0
#> 4    -1850 H3K27ac 90.92000 201.2272 16.43013    8.5
#> 5    -1800 H3K27ac 88.02667 192.6586 15.73051   10.0
#> 6    -1750 H3K27ac 85.52000 173.5301 14.16867   14.5
## we could then plot for instance using ggplot:
# ggplot(d, aes(position, mean, colour=sample)) + geom_line(size=1.2)
```
