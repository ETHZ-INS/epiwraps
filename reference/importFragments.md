# importFragments

A fast, lightweight replacement for `rtracklayer::import(format="bed")`
specialized for tabix-indexed fragment files (chr, start, end,
barcode\[, count\]).

## Usage

``` r
importFragments(x, which = NULL, useCLI = nzchar(Sys.which("tabix")))
```

## Arguments

- x:

  The path to a tabix-indexed fragment file, or a TabixFile object.

- which:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of the region(s) for which to read fragments. If missing, the whole
  file is read.

- useCLI:

  Logical; whether to stream the data through the \`tabix\` binary (much
  faster, if available) rather than \`rtracklayer::import\`.

## Value

A [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
with the same conventions as `rtracklayer::import` for BED files
(1-based, closed coordinates; unstranded), with the barcode in the
\`name\` column and, if present in the file, the read count in the
\`score\` column.

## Examples

``` r
frags <- tempfile(fileext = ".tsv")
d <- data.frame(chr=rep(letters[1:2], each=10), start=rep(100*(1:10),2))
d$end <- d$start + 15L
d$cell <- paste0("barcode",sample.int(3, nrow(d), replace=TRUE))
write.table(d, frags, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
frags <- Rsamtools::bgzip(frags)
Rsamtools::indexTabix(frags, format = "bed")
#> [1] "/tmp/RtmpCArGKG/file31321786f3ad.tsv.bgz.tbi"
importFragments(frags, GRanges("a", IRanges(1, 1000)))
#> GRanges object with 9 ranges and 1 metadata column:
#>       seqnames    ranges strand |        name
#>          <Rle> <IRanges>  <Rle> | <character>
#>   [1]        a   101-115      * |    barcode2
#>   [2]        a   201-215      * |    barcode3
#>   [3]        a   301-315      * |    barcode2
#>   [4]        a   401-415      * |    barcode2
#>   [5]        a   501-515      * |    barcode1
#>   [6]        a   601-615      * |    barcode3
#>   [7]        a   701-715      * |    barcode1
#>   [8]        a   801-815      * |    barcode2
#>   [9]        a   901-915      * |    barcode3
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
