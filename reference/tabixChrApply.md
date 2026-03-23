# tabixChrApply

Runs a function on reads/fragments from each chromosomes of a
Tabix-indexed fragment file. This is especially used by other functions
to avoid loading all alignments into memory, or to parallelize reads
processing.

## Usage

``` r
tabixChrApply(
  x,
  fn,
  keepSeqLvls = NULL,
  exclude = NULL,
  only = NULL,
  BPPARAM = NULL,
  progress = TRUE,
  ...
)
```

## Arguments

- x:

  The path to a tabix-indexed bam file, or a TabixFile object.

- fn:

  The function to be run, the first argument of which should be a
  \`GRanges\`

- keepSeqLvls:

  An optional vector of seqLevels to keep

- exclude:

  An optional GRanges of regions for which overlapping reads should be
  excluded.

- only:

  An optional GRanges of regions for which overlapping reads should be
  included. If set, all other reads are discarded.

- BPPARAM:

  A \`BiocParallel\` parameter object for multithreading. Note that if
  used, memory usage will be high; in this context we recommend a high
  \`nChunks\`.

- progress:

  Logical; whether to show a progress bar.

- ...:

  Passed to \`fn\`

## Value

A list of whatever \`fn\` returns

## Examples

``` r
# generate dummy regions and save them to a temp file:
frags <- tempfile(fileext = ".tsv")
d <- data.frame(chr=rep(letters[1:2], each=10), start=rep(100*(1:10),2))
d$end <- d$start + 15L
write.table(d, frags, col.names=FALSE, row.names=FALSE, sep="\t")
# tabix-index it
frags <- Rsamtools::bgzip(frags)
Rsamtools::indexTabix(frags, format = "bed")
#> [1] "/tmp/RtmpTktJrD/file389b267a9e40.tsv.bgz.tbi"
# now we can do something chunk-wise, e.g. extract coverage:
res <- tabixChrApply(frags, fn=coverage)
# aggregate the chunk results into an RleList object:
reduceRleLists(res)
#> RleList of length 2
#> $`"a"`
#> integer-Rle of length 1015 with 20 runs
#>   Lengths: 100  15  85  15  85  15  85  15 ...  85  15  85  15  85  15  85  15
#>   Values :   0   1   0   1   0   1   0   1 ...   0   1   0   1   0   1   0   1
#> 
#> $`"b"`
#> integer-Rle of length 1015 with 20 runs
#>   Lengths: 100  15  85  15  85  15  85  15 ...  85  15  85  15  85  15  85  15
#>   Values :   0   1   0   1   0   1   0   1 ...   0   1   0   1   0   1   0   1
#> 
```
