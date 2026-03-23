# peakPbCountsSE

Generate a pseudobulk peak counts SummarizedExperiment from a fragment
file.

## Usage

``` r
peakPbCountsSE(
  fragfile,
  peaks,
  bcmap,
  insertions = FALSE,
  genome = NULL,
  minFragLength = 1L,
  maxFragLength = 5000L
)
```

## Arguments

- fragfile:

  The path to a Tabix-indexed fragment file.

- peaks:

  A GRanges of the regions in which to count.

- bcmap:

  A named vector, indicating the pseudobulk sample (values) in which to
  include each barcode (names).

- insertions:

  If TRUE, (shifted) Tn5 insertions events are counted instead of
  fragments. This means that each fragment gets counted twice (for both
  ends). Default FALSE.

- genome:

  A optional genome object or path to a genom fasta file. If included,
  GC bias will be added to the rowData of the output object.

- minFragLength:

  Minimum fragment length for a fragment to be counted.

- maxFragLength:

  Maximum fragment length for a fragment to be counted.

## Value

A
[RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
with a 'counts' assay, and columns corresponding to each unique value of
\`bcmap\`.
