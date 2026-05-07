# frag2bw

Creates a coverage bigwig file from a Tabix-indexed fragment file.

## Usage

``` r
frag2bw(
  tabixFile,
  output_bw,
  binWidth = 20L,
  scaling = TRUE,
  type = c("full", "center", "start", "end", "ends"),
  barcodes = NULL,
  strand = c("*", "+", "-"),
  shift = 0L,
  log1p = FALSE,
  exclude = NULL,
  minFragLength = 1L,
  maxFragLength = 5000L,
  keepSeqLvls = NULL,
  useScore = FALSE,
  forceSeqlevelsStyle = NULL,
  only = NULL,
  format = "bed",
  binSummarization = c("mean", "max", "min"),
  verbose = TRUE
)
```

## Arguments

- tabixFile:

  The path to a tabix-indexed bam file, or a TabixFile object.

- output_bw:

  The path to the output bigwig file

- binWidth:

  The window size. A lower value (min 1) means a higher resolution, but
  larger file size.

- scaling:

  Either TRUE (performs Count Per Million scaling), FALSE (no scaling),
  or a numeric value by which the signal will be divided. If \`bgbam\`
  is given and \`scaling=TRUE\`, the background will be scaled to the
  main signal.

- type:

  Type of the coverage to compile. Either full (full read/fragment),
  start (count read/fragment start locations), end, center, or 'ends'
  (both ends of the read/fragment).

- barcodes:

  An optional list of barcodes to use (assuming that the file contains
  the column)

- strand:

  Strand(s) to capture (any by default).

- shift:

  Shift (from 3' to 5') by which reads/fragments will be shifted. If
  \`shift\` is an integer vector of length 2, the first value will
  represent the shift for the positive strand, and the second for the
  negative strand.

- log1p:

  Whether to log-transform (\`log(x+1)\`) the (scaled) signal.

- exclude:

  An optional GRanges of regions for which overlapping reads should be
  excluded.

- minFragLength:

  Minimum fragment length (ignored if \`paired=FALSE\`)

- maxFragLength:

  Maximum fragment length (ignored if \`paired=FALSE\`)

- keepSeqLvls:

  An optional vector of seqLevels (i.e. chromosomes) to include.

- useScore:

  Whether to use the score column (if any) as coverage weights.

- forceSeqlevelsStyle:

  If specified, forces the use of the specified seqlevel style for the
  output bigwig. Can take any value accepted by \`seqlevelsStyle\`.

- only:

  An optional GRanges of regions for which overlapping reads should be
  included. If set, all other reads are discarded.

- format:

  The format of the fragment file.

- binSummarization:

  The method to summarize nucleotides into each bin, either "mean"
  (default), "min" or "max".

- verbose:

  Logical; whether to print progress messages

## Value

The bigwig filepath. Alternatively, if \`output_bw=NA_character\_\`, the
coverage data is not written to file but returned.

## Examples

``` r
# we first create a fake tabix file:
library(GenomicRanges)
library(rtracklayer)
reads <- GRanges(rep(c("1","2"), c(5,2)),
                 IRanges(5000+10*1:7, width=100))
bedf <- tempfile(fileext=".bed")
rtracklayer::export.bed(reads, bedf)
bedf <- Rsamtools::bgzip(bedf)
Rsamtools::indexTabix(bedf, format="bed")
#> [1] "/tmp/RtmpG77iEa/file1c2643a83f5.bed.bgz.tbi"
# convert to bigwig
frag2bw(bedf, tempfile(fileext=".bw"))
#> Reading in signal...
#> Writing bigwig...
```
