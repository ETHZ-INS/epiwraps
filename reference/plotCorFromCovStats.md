# plotCorFromCovStats

plotCorFromCovStats

## Usage

``` r
plotCorFromCovStats(
  qc,
  method = c("pearson", "spearman"),
  col = NULL,
  na_col = "white",
  column_title = NULL,
  ...
)
```

## Arguments

- qc:

  A list of coverage statistics, as produced by
  [`getCovStats`](https://ethz-ins.github.io/epiwraps/reference/getCovStats.md).

- method:

  The correlation metrics to include

- col:

  Optional heatmap colors

- na_col:

  Color for the diagonal; passed to
  [`Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

- column_title:

  Column title (if NULL, uses the metric)

- ...:

  Passed to
  [`Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

## Value

A \`Heatmap\` or \`HeatmapList\` object ready to be plotted.
