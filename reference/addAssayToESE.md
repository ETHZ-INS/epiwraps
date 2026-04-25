# addAssayToESE

Adds an assay of signal matrices to an existing \`EnrichmentSE\` object.

## Usage

``` r
addAssayToESE(x, a, name = "normalized", replace = TRUE)
```

## Arguments

- x:

  An object of class \`EnrichmentSE\`, as produced by
  [`signal2Matrix`](https://ethz-ins.github.io/epiwraps/reference/signal2Matrix.md).

- a:

  The assay to add, e.g. a list of normalizedMatrix objects

- name:

- replace:

  Logical, whether to replace any existing assay of the same name
  (default TRUE). If FALSE and the assay already existed, the new assay
  name is given a suffix.

## Value

\`x\` with the added/updated assay.

## Examples

``` r
# we first get an EnrichmentSE object:
data(exampleESE)
# then we will create a new assay which is simply sqrt-transformed, and add 
# it back in the object
newAssay <- lapply(getSignalMatrices(exampleESE), sqrt)
exampleESE <- addAssayToESE(exampleESE, newAssay, name="sqrt")
```
