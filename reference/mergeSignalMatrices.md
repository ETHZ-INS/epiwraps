# mergeSignalMatrices: aggregates two or more signal matrices.

mergeSignalMatrices: aggregates two or more signal matrices.

## Usage

``` r
mergeSignalMatrices(ml, aggregation = c("mean", "sum", "median"), assay = 1L)
```

## Arguments

- ml:

  A named list of signal matrices or an EnrichmentSE object as produced
  by
  [`signal2Matrix`](https://ethz-ins.github.io/epiwraps/reference/signal2Matrix.md)

- aggregation:

  The method to aggregate matrices

- assay:

  Assay to use (ignored unless \`ml\` is an ESE object), defaults to the
  first assay.

## Value

A single \`normalizedMatrix\` object.

## Examples

``` r
# we first get an EnrichmentSE object:
data(exampleESE)
# we merge the two tracks:
merged <- mergeSignalMatrices(exampleESE)
# we could then plot the merge (not run):
# plotEnrichedHeatmaps(merged)
```
