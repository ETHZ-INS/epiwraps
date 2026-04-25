# estimateFragSize

estimateFragSize

## Usage

``` r
estimateFragSize(
  bam,
  ctrl = NULL,
  binSize = 10L,
  mfold = c(10, 50),
  ...,
  minSummitCount = 8L,
  useSeqLevels = NULL,
  maxSize = 2500L,
  priorLength = 200L,
  blacklist = NULL,
  ret = c("mode", "median", "mean", "tables", "plots", "distances"),
  BPPARAM = SerialParam()
)
```

## Arguments

- bam:

  The path to one or more bam files

- ctrl:

  Optional path to a control bam file (if \`length(bam)\>1\`, the same
  control will be used for all).

- binSize:

  Bin size. The precision of the reported fragment size is necessary
  lower than this. A higher bin size will improve the summit
  identification in low-coverage regions. We recommend leaving the
  default value.

- mfold:

  The range of fold-enrichment over the control (if \`ctrl\` provided)
  or of coverages for the identification of regions based on which
  distance will be estimated.

- minSummitCount:

  The minimum read count for a summit to be considered.

- useSeqLevels:

  An optional vector of seqLevels in which to conduct the analysis.

- maxSize:

  The maximum size of regions to be used

- priorLength:

  The prior fragment length (use for read extension to identify enriched
  regions)

- blacklist:

  Optional \`GRanges\` of blacklisted regions to be excluded.

- ret:

  The type of return, either a 'table' of pairs of summits and their
  properties, or a 'plot', or the median/mean/mode of the distance
  distribution.

- BPPARAM:

  A \`BiocParallel\` parameter object for multithreading. Only used if
  multiple files are given in \`bam\`.

## Value

By default, the estimated (mode) fragment length(s), but see the \`ret\`
argument

## Examples

``` r
# get an example bam file
bam <- system.file("extdata", "ex1.bam", package="Rsamtools")
suppressWarnings(estimateFragSize(bam))
#> [1] 1433
```
