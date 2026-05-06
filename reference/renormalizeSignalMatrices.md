# renormalizeSignalMatrices

Renormalizes a list of signal matrices or an EnrichmentSE object.

## Usage

``` r
renormalizeBorders(ml, trim = NULL, assay = "input", nWindows = NULL)

renormalizeSignalMatrices(
  ml,
  method = c("border", "top", "manual"),
  trim = NULL,
  fromAssay = "input",
  toAssay = NULL,
  nWindows = NULL,
  scaleFactors = NULL
)
```

## Arguments

- ml:

  A named matrix list or EnrichmentSE object as produced by
  [`signal2Matrix`](https://ethz-ins.github.io/epiwraps/reference/signal2Matrix.md).

- trim:

  Quantiles trimmed at each extreme before calculating normalization
  factors.

- assay:

  The name of the assay to use as input.

- nWindows:

  Number of border windows/bins to use for border normalization.

- method:

  Either "border" or "top" (see details below).

- fromAssay:

  Assay to use (ignored unless \`ml\` is an EnrichmentSE object),
  defaults to the first assay.

- toAssay:

  Assay in which to store the normalized data (ignored unless \`ml\` is
  an EnrichmentSE object). By default an assay name will be set based on
  the normalization method used.

- scaleFactors:

  A numeric vector of same length as \`ml\`, indicating the scaling
  factors by which to multiply each matrix. Alternatively, a numeric
  matrix with a number of rows equal to the length of \`ml\`, and two
  columns indicating the alpha and beta arguments of a s3norm
  normalization. Ignored unless \`method="manual"\`.

## Value

Either a renormalized list of signal matrices or, if \`ml\` was an
\`EnrichmentSE\` object, the same object with an additional normalized
assay automatically put at the front.

## Details

\* \`method="border"\` works on the assumption that the left/right
borders of the matrices represent background signal which should be
equal across samples. As a result, it will work only if 1) the
left/right borders of the matrices are sufficiently far from the signal
(e.g. peaks) to be chiefly noise, and 2) the signal-to-noise ratio is
comparable across tracks/samples. \* \`method="top"\` instead works on
the assumption that the highest signal should be the same across
tracks/samples. By default, extreme values are trimmed before
establishing either kind of normalization factor. The proportion trimmed
can be set using the \`trim\` argument, and is by default 10 \*
\`method="manual"\` enables the use of independently computed
normalization factors, for instance obtained through
[`getNormFactors`](https://ethz-ins.github.io/epiwraps/reference/getNormFactors.md).

## Functions

- `renormalizeBorders()`: deprecated \> renormalizeSignalMatrices

## Examples

``` r
# we first get an EnrichmentSE object:
data(exampleESE)
# we normalize them
exampleESE <- renormalizeSignalMatrices(exampleESE)
# see the `vignette("multiRegionPlot")` for more info on normalization.
```
