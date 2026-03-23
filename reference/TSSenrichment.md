# TSSenrichment

A common quality metric for ATAC-seq data (see
\[definition\](https://www.encodeproject.org/data-standards/terms/#enrichment)
). The interpretation of the enrichment score depends on the annotation
(the score tends to increase when only fewer, more common TSS are used),
but according to ENCODE guidelines anything below 5 is of concerning
quality, while a score \>8 is ideal.

## Usage

``` r
TSSenrichment(tracks, ensdb, useSeqLevels = NULL)
```

## Arguments

- tracks:

  A (named) vector of paths to bigwig files.

- ensdb:

  An \`ensembldb\` object. Alternatively, a GRanges object of regions
  centered around TSS.

- useSeqLevels:

  Optional seqlevels to use. If NULL, all are used.

## Value

A list with the slots \`score\` (numeric vector of TSS enrichment scores
per sample) and \`data\` (per bin enrichment, for plotting)

## Examples

``` r
# we first fetch the path to the example bigwig file:
bw <- system.file("extdata/example_atac.bw", package="epiwraps")
## normally, we would load an ensembldb object using AnnotationHub. For the 
## purpose of this example, we'll pretend that the following set of regions
## represent TSS:
tss <- system.file("extdata/example_peaks.bed", package="epiwraps")
tss <- rtracklayer::import(tss)
en <- TSSenrichment(bw, tss)
#> Reading /home/runner/work/_temp/Library/epiwraps/extdata/example_atac.bw
en$score
#> example_atac 
#>     5.768383 
## you can also plot using something like this:
## ggplot(en$data, aes(position, enrichment, colour=sample)) + geom_line()
```
