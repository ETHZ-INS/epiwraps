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
# we'll make a lightweight ensembldb annotation for the annotation:
library(ensembldb)
#> Loading required package: GenomicFeatures
#> Loading required package: AnnotationDbi
#> Loading required package: AnnotationFilter
#> 
#> Attaching package: 'ensembldb'
#> The following object is masked from 'package:stats':
#> 
#>     filter
chrY <- system.file("chrY", package="ensembldb")
edb <- EnsDb(makeEnsemblSQLiteFromTables(path=chrY ,dbname=tempfile()))
#> Processing 'chromosome' table ... 
#> OK
#> Processing 'gene' table ... 
#> OK
#> Processing 'trancript' table ... 
#> OK
#> Processing 'exon' table ... 
#> OK
#> Processing 'tx2exon' table ... 
#> OK
#> Creating indices ... 
#> OK
#> Checking validity of the database ... 
#> OK
# we run teh annotation:
regions <- annotateRegions(regions, edb)
# this adds metadata columns to the regions:
regions
#> GRanges object with 2 ranges and 5 metadata columns:
#>       seqnames            ranges strand | distance2nearestTSS
#>          <Rle>         <IRanges>  <Rle> |           <integer>
#>   [1]     chrY   2655742-2655890      * |                   2
#>   [2]     chrY 28730110-28730950      * |               -6798
#>       nearestTSS.gene_name      nearestTSS TSS.overlap             class
#>                <character>     <character>    <factor>          <factor>
#>   [1]                  SRY ENST00000383070  intergenic proximal <=1000bp
#>   [2]           SLC25A15P1 ENST00000456738  intergenic intergenic       
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
table(regions$class)
#> 
#>        intergenic proximal <=1000bp 
#>                 1                 1 
```
