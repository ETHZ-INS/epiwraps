# annotateRegions

Annotates a GRanges on the basis of an annotation object (e.g.
[`EnsDb`](https://rdrr.io/pkg/ensembldb/man/EnsDb.html)).

## Usage

``` r
annotateRegions(
  regions,
  anno,
  proximal = c(2500, 1000),
  filter = AnnotationFilterList(),
  extra = list(),
  ignore.strand = TRUE,
  ...
)
```

## Arguments

- regions:

  A GRanges object

- anno:

  An annotation object, such as an
  [`EnsDb`](https://rdrr.io/pkg/ensembldb/man/EnsDb.html)) object, a
  TxDb object, or a GRanges object of a GENCODE-like gtf or the path to
  such a file.

- proximal:

  The threshold(s) for TSS proximal regions. Multiple values will result
  in multiple class factor levels.

- filter:

  An
  [`AnnotationFilter`](https://rdrr.io/pkg/AnnotationFilter/man/AnnotationFilter.html)
  to filter transcripts. Only used if \`anno\` is an
  [`EnsDb`](https://rdrr.io/pkg/ensembldb/man/EnsDb.html)).

- extra:

  An optional named list of GRanges for additional overlaps. Each list
  element will create an additional binary metadata column.

- ignore.strand:

  Whether to ignore the strand for the overlap with elements of
  \`extra\` (default TRUE).

- ...:

  Passed to
  [`overlapsAny`](https://rdrr.io/pkg/GenomicRanges/man/findOverlaps-methods.html)
  for the overlaps with \`extra\`.

## Value

The sorted \`regions\` object with additional annotation columns.

## Examples

``` r
# we first create some regions we want to annotate:
regions <- as(c("chrY:2655742-2655890", "chrY:28730110-28730950"), "GRanges")
# we'll use a lightweight ensembldb annotation for the annotation:
chrY <- system.file("chrY", package="ensembldb")
edb <- EnsDb(makeEnsemblSQLiteFromTables(path=chrY ,dbname=tempfile()))
#> Error in EnsDb(makeEnsemblSQLiteFromTables(path = chrY, dbname = tempfile())): could not find function "EnsDb"
# we run teh annotation:
regions <- annotateRegions(regions, edb)
#> Error: object 'edb' not found
# this adds metadata columns to the regions:
regions
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames            ranges strand
#>          <Rle>         <IRanges>  <Rle>
#>   [1]     chrY   2655742-2655890      *
#>   [2]     chrY 28730110-28730950      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
table(regions$class)
#> < table of extent 0 >
```
