# motifFootprint

This is a wrapper around \`ATACseqQC::factorFootprints\`, making it
compatible with a broader variety of inputs and use cases.

## Usage

``` r
motifFootprint(bamfile, motif, motif_occurences, genome = NULL, around = 100)
```

## Arguments

- bamfile:

  The path to a bam file

- motif:

  The motif (as probability matrix)

- motif_occurences:

  A GRanges of the motif occurences

- genome:

  Optional genome used for seqlengths (otherwise estimated from the
  occurences)

- around:

  How much around the motif to plot

## Value

A plot
