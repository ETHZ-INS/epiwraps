# peakCountsFromBAM

Creates a SummarizedExperiment of fragment (or insertion) counts from
bam files that overlap given regions.

## Usage

``` r
peakCountsFromBAM(
  bam_files,
  regions,
  paired,
  extend = 0L,
  shift = 0L,
  type = c("full", "center", "start", "end", "ends"),
  ov.type = "any",
  maxgap = -1L,
  minoverlap = 1L,
  ignore.strand = TRUE,
  strandMode = 1,
  includeDuplicates = TRUE,
  includeSecondary = FALSE,
  minMapq = 1L,
  minFragLength = 1L,
  maxFragLength = 5000L,
  splitByChr = 3,
  randomAcc = FALSE,
  getMedianFragLength = FALSE,
  verbose = TRUE
)
```

## Arguments

- bam_files:

  A vector of paths to the bam files.

- regions:

  A \`GRanges\` of regions in which to counts.

- paired:

  Logical; whether the data is paired (assumed unpaired by default). Use
  \`paired="auto"\` for automatic detection using the first bam file.

- extend:

  The amount \*by\* which to extend single-end reads (e.g. fragment
  length minus read length). If \`paired=TRUE\` and \`type\` is either
  'ends' or 'center', then the extension will be applied after taking
  the (shifted) fragment ends or centers, resulting in ranges of width
  equal to \`extend\`.

- shift:

  Shift (from 3' to 5') by which reads/fragments will be shifted. If
  \`shift\` is an integer vector of length 2, the first value will
  represent the shift for the positive strand, and the second for the
  negative strand.

- type:

  Type of the coverage to compile. Either full (full read/fragment),
  start (count read/fragment start locations), end, center, or 'ends'
  (both ends of the read/fragment).

- ov.type:

  Overlap type. See the \`type\` argument of
  `link[GenomicRanges]{countOverlaps}`.

- maxgap:

  Maximum gap allowed for overlaps (see the corresponding argument of
  `link[GenomicRanges]{countOverlaps}`).

- minoverlap:

  Minimum overlap (see the corresponding argument of
  `link[GenomicRanges]{countOverlaps}`).

- ignore.strand:

  Logical; whether to ignore strand for the purpose of counting overlaps
  (default TRUE).

- strandMode:

  The strandMode of the data (whether the strand is given by the first
  or second mate, which depends on the library prep protocol). See
  [strandMode](https://rdrr.io/pkg/GenomicAlignments/man/GAlignmentPairs-class.html)
  for more information. This parameter has no effect unless one of the
  \`strand\`, \`extend\` parameters or a strand-specific \`shift\` are
  used.

- includeDuplicates:

  Logical, whether to include reads flagged as duplicates.

- includeSecondary:

  Logical; whether to include secondary alignments

- minMapq:

  Minimum mapping quality (1 to 255)

- minFragLength:

  Minimum fragment length (ignored if \`paired=FALSE\`)

- maxFragLength:

  Maximum fragment length (ignored if \`paired=FALSE\`)

- splitByChr:

  Whether to process chromosomes separately, and if so by how many
  chunks. The should not affect the output, and is simply slightly
  slower and consumes less memory. Can be a logical value (in which case
  each chromosome is processed separately), but we instead recommend
  giving a positive integer indicating the number of chunks.

- randomAcc:

  Logical, whether to use random access. This is disabled by default
  because the overhead of random access to a lot of regions is typically
  worse than reading the entire file. However, if you need to get counts
  in few regions, enabling this will be faster. Note however that when
  using random access, the output object will not contain depth
  information.

- getMedianFragLength:

  Logical; whether to compile the median fragment length per region.
  This is slightly slower. The log10-transformed, (weighted mean across
  samples of the) median fragment length per region is stored in
  \`rowData(results)\$flbias\`.

- verbose:

  Logical; whether to print progress messages

## Value

A
[`RangedSummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
with a 'counts' assay.

## Examples

``` r
# get an example bam file
bam <- system.file("extdata", "ex1.bam", package="Rsamtools")
# create regions of interest
peaks <- GRanges(c("seq1","seq1","seq2"), IRanges(c(400,900,500), width=100))
peakCountsFromBAM(bam, peaks, paired=FALSE)
#> Reading file /home/runner/work/_temp/Library/Rsamtools/extdata/ex1.bam
#> class: RangedSummarizedExperiment 
#> dim: 3 1 
#> metadata(0):
#> assays(1): counts
#> rownames(3): seq1:400-499 seq1:900-999 seq2:500-599
#> rowData names(0):
#> colnames(1): ex1
#> colData names(2): total_depth depth
```
