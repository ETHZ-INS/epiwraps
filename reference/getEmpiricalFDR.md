# getEmpiricalFDR

Computes the FDR as the proportion of negative peaks with a more extreme
p-value, eventually extrapolating and preserving the ranking.

## Usage

``` r
getEmpiricalFDR(log10p, pneg, n = length(log10p) * 10)
```

## Arguments

- log10p:

  -log10 p-values of candidate peaks

- pneg:

  -log10 p-values of negative peaks

- n:

  The number of hypotheses (used when the empirical FDR is zero)

## Value

A data.frame with the empirical FDR and a smoothed -log10(FDR)
