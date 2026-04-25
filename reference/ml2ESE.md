# Creates an EnrichmentSE from a list of normalizedMatrix objects

Creates an EnrichmentSE from a list of normalizedMatrix objects

## Usage

``` r
ml2ESE(ml, rowRanges, assayName = "input", addScore = FALSE, ...)
```

## Arguments

- ml:

  A named list of normalizedMatrix objects with corresponding rows.

- rowRanges:

  An optional GRanges object corresponding to the rows of each object of
  \`ml\`.

- assayName:

  The name of the assay, defaults to 'input'

- addScore:

  Logical; whether to add an enriched_score assay.

- ...:

  Passed to \`SummarizedExperiment()\`

## Value

An \`EnrichedSE\` object, inheriting from a RangedSummarizedExperiment.

## Examples

``` r
# for an example we first need a list of signal matrices. To this end, 
# we first fetch the path to the example bigwig file:
bw <- system.file("extdata/example_atac.bw", package="epiwraps")
# we load example regions:
regions <- rtracklayer::import(system.file("extdata/example_peaks.bed", 
                                           package="epiwraps"))
# we obtain the matrix of the signal around the regions, indicating that we
# want the output as a list of signal matrices:
m <- signal2Matrix(bw, regions, ret="list")
#> Reading /home/runner/work/_temp/Library/epiwraps/extdata/example_atac.bw
# we can then transform this into an EnrichmentSE object:
m <- ml2ESE(m, rowRanges=regions)
```
