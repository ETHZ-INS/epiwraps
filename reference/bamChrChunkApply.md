# bamChrChunkApply

Runs a function on reads/fragments from chunks of (chromosomes) of an
indexed bam file. This is especially used by other functions to avoid
loading all alignments into memory, or to parallelize reads processing.

## Usage

``` r
bamChrChunkApply(
  x,
  FUN,
  paired = FALSE,
  keepSeqLvls = NULL,
  nChunks = 4,
  strandMode = 2,
  flgs = scanBamFlag(),
  exclude = NULL,
  mapqFilter = NA_integer_,
  progress = TRUE,
  BPPARAM = NULL,
  ...
)
```

## Arguments

- x:

  A bam file.

- FUN:

  The function to be run, the first argument of which should be a
  \`GRanges\`

- paired:

  Logical; whether to consider the reads as paired (fragments, rather
  than reads, will be returned)

- keepSeqLvls:

  An optional vector of seqLevels to keep

- nChunks:

  The number of chunks to use (higher will use less memory but increase
  overhead)

- strandMode:

  Strandmode for paired data (see
  [`readGAlignmentPairs`](https://rdrr.io/pkg/GenomicAlignments/man/readGAlignments.html)).

- flgs:

  \`scanBamFlag\` for filtering the reads

- exclude:

  An optional GRanges of regions for which overlapping reads should be
  excluded.

- mapqFilter:

  Integer of the minimum mapping quality for reads to be included.

- progress:

  Logical; whether to show a progress bar.

- BPPARAM:

  A \`BiocParallel\` parameter object for multithreading. Note that if
  used, memory usage will be high; in this context we recommend a high
  \`nChunks\`.

- ...:

  Passed to \`FUN\`

## Value

A list of whatever \`FUN\` returns

## Examples

``` r
# as an example we'll use the function to obtain fragment sizes:
bam <- system.file("extdata", "ex1.bam", package="Rsamtools")
fragLen <- bamChrChunkApply(bam, paired=TRUE, FUN=function(x) width(x))
quantile(unlist(fragLen))
#>   0%  25%  50%  75% 100% 
#>   54  199  209  218  243 
```
