# exportNarrowPeaks

exportNarrowPeaks

## Usage

``` r
exportNarrowPeaks(p, file)
```

## Arguments

- p:

  A GRanges, as produced by \`callPeaks\`

- file:

  The path to the file where to save

## Value

Invisible file path.

## Examples

``` r
# call some peaks:
bam <- system.file("extdata", "ex1.bam", package="Rsamtools")
peaks <- callPeaks(bam, paired=TRUE)
#> Reading signal and identifying candidate regions...
#> Identified 2 candidate regions
#> Computing significance...
#> (In the absence of a control, FDR is unlikely to be calibrated)
#> Reporting 2 regions, 2 with FDR<0.05
# save them:
filepath <- tempfile(fileext="narrowPeak")
exportNarrowPeaks(peaks, filepath)
```
