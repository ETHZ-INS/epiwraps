# signalsAcrossSamples

Obtain a matrix of score/coverages across a region for a list of BigWig
files.

## Usage

``` r
signalsAcrossSamples(files, region, ignore.strand = TRUE)
```

## Arguments

- files:

  A named list of paths to biwgig files or of \`GRanges\` objects with a
  \`score\` column.

- region:

  The region of interest, either given as a string (in the
  "chr:start-end" format) or as a \`GRanges\` of length 1.

- ignore.strand:

  Logical; whether to merge scores from the two strands given stranded
  objects.

## Value

A disjoined \`GRanges object\` with the scores as metadata columns.
