# peakCountsFromFrags

Creates a SummarizedExperiment of cell-level or pseudo-bulk-level
fragment (or insertion) counts from a tabix-indexed fragment file and a
set of regions of interest.

## Usage

``` r
peakCountsFromFrags(
  fragfile,
  regions,
  barcodes = NULL,
  insertions = FALSE,
  minFragLength = 1L,
  maxFragLength = 5000L,
  ov.type = "any",
  maxgap = -1L,
  minoverlap = 1L,
  ignore.strand = TRUE,
  ...
)
```

## Arguments

- fragfile:

  A path to the tabix-indexed fragment file.

- regions:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of regions in which to count overlaps.

- barcodes:

  An optional character vector of cell barcodes to include. If provided,
  only these barcodes will be considered, if \`NULL\`, all barcodes
  included in the file are used. If \`barcodes\` is a named vector, the
  names will be considered to represent the cell barcode, the values to
  represent the pseudo-bulk sample in which to include the respective
  barcodes, and pseudo-bulk counts will be returned.

- insertions:

  Logical; if `TRUE`, the ends of the fragments (insertions) are counted
  instead of the entire fragment. This means each fragment can
  contribute up to two counts. Default `FALSE`.

- minFragLength:

  Minimum fragment length to be considered. Default 1.

- maxFragLength:

  Maximum fragment length to be considered. Default 5000.

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

- ...:

  Passed to
  [`tabixChrApply`](https://ethz-ins.github.io/epiwraps/reference/tabixChrApply.md).

## Value

A
[`RangedSummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
with a 'counts' assay. The mean fragment length per region is also
stored in the \`rowData\` of the object.

## Examples

``` r
# generate dummy regions and save them to a temp file:
frags <- tempfile(fileext = ".tsv")
d <- data.frame(chr=rep(letters[1:2], each=10), start=rep(100*(1:10),2))
d$end <- d$start + 15L
d$cell <- paste0("barcode",sample.int(3, nrow(d), replace=TRUE))
write.table(d, frags, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
# tabix-index it
frags <- Rsamtools::bgzip(frags)
Rsamtools::indexTabix(frags, format = "bed")
#> [1] "/tmp/Rtmpwo4ZtI/file1c465b4788d0.tsv.bgz.tbi"
# we create regions of interest:
regions <- GRanges(c("a","b"), IRanges(400,width=300))
# we get the counts:
se <- peakCountsFromFrags(frags, regions)
se
#> class: RangedSummarizedExperiment 
#> dim: 2 3 
#> metadata(0):
#> assays(1): counts
#> rownames(2): a:400-699 b:400-699
#> rowData names(1): flbias
#> colnames(3): barcode1 barcode2 barcode3
#> colData names(0):
# we could also get pseudobulk counts by passing a barcode map:
bcmap <- setNames(c("PB1","PB1","PB2"),paste0("barcode",1:3))
pb <- peakCountsFromFrags(frags, regions, barcodes=bcmap)
pb
#> class: RangedSummarizedExperiment 
#> dim: 2 2 
#> metadata(0):
#> assays(1): counts
#> rownames(2): a:400-699 b:400-699
#> rowData names(1): flbias
#> colnames(2): PB1 PB2
#> colData names(1): depth
```
