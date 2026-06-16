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
if (FALSE) { # \dontrun{
# assuming we have a bigwig file and an ensdb object:
en <- TSSenrichment(bw, ensb)
en$score
# you can also plot using something like this:
ggplot(en$data, aes(position, enrichment, colour=sample)) + geom_line()
} # }
```
