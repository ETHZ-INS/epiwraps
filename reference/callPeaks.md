# callPeaks

This is a native R peak caller loosely based on the general MACS2/3
strategy (Zhang et al., Genome Biology 2008). It can work from BAM
files, fragment files, or from coverage tracks. The results are highly
concordant with MACS, although the significance estimates are slightly
more conservative, the peak resolution differs, and the function
requires much more memory (see details). We believe its calls are
slightly better than MACS3 in the absence of a control, and slightly
worse when using a control.

## Usage

``` r
callPeaks(
  bam,
  ctrl = NULL,
  paired,
  type = c("narrow", "broad"),
  fragLength = 200L,
  globalNullH = FALSE,
  gsize = NULL,
  blacklist = NULL,
  binSize = 10L,
  nChunks = 8L,
  flags = scanBamFlag(isDuplicate = FALSE),
  minPeakCount = 5L,
  minFoldEnr = 1.4,
  pthres = 10^-5,
  maxSize = NULL,
  bgWindow = c(1, 5, 10) * 1000,
  pseudoCount = 0.5,
  useStrand = !paired,
  outFormat = c("narrowPeak", "custom"),
  verbose = TRUE,
  ...
)
```

## Arguments

- bam:

  A signal bam file (which should be accompanied by an index file).
  Alternatively, a TABIX-indexed fragment file, or an RleList object.

- ctrl:

  An optional path to a control bam file. Alternatively, an RleList
  object.

- paired:

  Logical, whether the reads are paired.

- type:

  The type of peaks to identify ('narrow' or 'broad').

- fragLength:

  Fragment length. Ignored if \`paired=TRUE\`. If \`useStrand=TRUE\`
  (default), this is only used for the initial candidate region
  identification, and sizes are adjusted after, so it doesn't need to be
  very precise.

- globalNullH:

  Logical; whether to use a global expectation, rather than the local
  background (default), in the absence of a control.

- gsize:

  The mappable genome size. Ignored unless \`globalNullH=TRUE\`. Can
  also be a species acronym in 'hs', 'mm', 'dm', and 'ce'.

- blacklist:

  An optional \`GRanges\` of regions to be excluded (or the path to such
  a file). Since the blacklisted regions are removed from both the
  signal and control peaks, this also has an important impact on the
  empirical FDR (when \`ctrl\` is given).

- binSize:

  Binsize used to estimate peak shift.

- nChunks:

  The number of processing chunks. Increasing this will reduce memory
  usage.

- flags:

  An optional
  [scanBamFlag](https://rdrr.io/pkg/Rsamtools/man/ScanBamParam-class.html)
  object to filter the reads. By default, reads flagged as optical
  duplicates are excluded.

- minPeakCount:

  The minimum summit count for a region to be considered. Decreasing
  this can substantially increase the running time.

- minFoldEnr:

  The minimum fold-enrichment for a region to be considered. Decreasing
  this can substantially increase the running time.

- pthres:

  The p-value threshold to use.

- maxSize:

  The loose maximum size of a peak. This is the size above which the
  method will attempt to break up peaks into smaller ones. The default
  is 5000 for \`type="broad"\`, and will depend on the estimated
  fragment size for narrow peak calls.

- bgWindow:

  The windows to consider (in addition to the peak itself) for local
  background.

- pseudoCount:

  The pseudocount to use when computing logFC.

- useStrand:

  Logical; whether to use strand information to better estimate the peak
  boundaries with single-end data.

- outFormat:

  The output format ('custom' or 'narrowPeak')

- verbose:

  Logical; whether to output progress messages

- ...:

  Passed to
  [`bamChrChunkApply`](https://ethz-ins.github.io/epiwraps/reference/bamChrChunkApply.md)

## Value

A \`GRanges\`, by default following the narrowPeak format specification.

## Details

Unless \`globalNullH=TRUE\`, this function uses MACS' local lambda
(defined by \`bgWindow\`). A major difference with MACS2/3 is that,
rather than using sliding windows, if works on the running list encoding
of the coverage(s). As a consequence, significance is estimated based on
the peak's maximum coverage, which is very similar for narrow peaks, but
very different for broad peaks, which will not produce astronomical
p-values as is the case with MACS. Peak sizes will also differ and the
distribution tends to be narrower (i.e. fewer small and very large
peaks).

Because FDR calculation is so tricky for peak calling, the function uses
an uncorrected p-value threshold for its output. If desired, users can
further filter based on the FDR. We recommend including a \`blacklist\`
for FDR estimation.

On a single thread, the function is nearly as fast as MACS2/3, but uses
much more memory, chiefly due to the fact that reads are read into
memory in large chunks. On a single thread, memory usage will be roughly
the sum of the filesize of your bam files (i.e. IP + input, if using an
input) divided by \`nChunks-1\`.

The function uses
[`bamChrChunkApply`](https://ethz-ins.github.io/epiwraps/reference/bamChrChunkApply.md)
to obtain the coverages, and can accept any argument of that function.
This means for instance that the \`mapqFilter\` argument can be used to
restrict the reads used, or a \`BPPARAM\` argument to enable
multi-threading (beware of memory consumption!).

To export as a narrowPeak file, see
[`exportNarrowPeaks`](https://ethz-ins.github.io/epiwraps/reference/exportNarrowPeaks.md).

## Examples

``` r
# we use the example bam file from the Rsamtools package:
bam <- system.file("extdata", "ex1.bam", package="Rsamtools")
peaks <- callPeaks(bam, paired=TRUE)
#> Reading signal and identifying candidate regions...
#> Identified 2 candidate regions
#> Computing significance...
#> (In the absence of a control, FDR is unlikely to be calibrated)
#> Reporting 0 regions, 0 with FDR<0.05
```
