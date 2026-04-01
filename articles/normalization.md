# Normalizing genomic signals

Abstract

This vignette covers the functions for normalizing genomic signals.
Since this is illustrated with visualization, it is recommended that you
read the `bam2bw` and `multiRegionPlot` vignettes first.

## Introduction

`epiwraps` includes two ways of calculating normalization factors:
either from the signal files (e.g. bam or bigwig files), which is the
most robust way and enables all options, or from an `EnrichmentSE`
object (see the [multiRegionPlot
vignette](https://ethz-ins.github.io/epiwraps/articles/multiRegionPlot.md)
for an intro to such an object) or signal matrices. In both cases, the
logic is the same: we estimate normalization factors (mostly single
linear scaling factors, although some methods involve more complex
normalization), and then apply them to signals that were extracted using
[`signal2Matrix()`](https://ethz-ins.github.io/epiwraps/reference/signal2Matrix.md).

### Applying normalization factors when generating the bigwig files

It is possible to also directly use computed normalization factors when
creating bigwig files. By default, the
[`bam2bw() function`](https://ethz-ins.github.io/epiwraps/articles/bam2bw.md)
scales using library size, which can be disabled using `scaling=FALSE`.
However, it is also possible to pass the `scaling` argument a manual
scaling factor, as computed by the functions described here. In this
vignette, however, we will focus on normalizing signal matrices.

## Obtaining normalization factors for a set of signal files

The
[`getNormFactors()`](https://ethz-ins.github.io/epiwraps/reference/getNormFactors.md)
function can be used to estimate normalization factors from either bam
or bigwig files. The files cannot be mixed (bam/bigwig), however, and it
is important to note that *normalization factors calculated on bam files
cannot be applied to data extracted from bigwig files, or vice versa*,
because the bigwig files are by default already normalized for library
size. If needed, however,
[`getNormFactors()`](https://ethz-ins.github.io/epiwraps/reference/getNormFactors.md)
can be used to apply the same method to both kind of files.

### Normalization methods

Simple library size normalization, as done by
[`bam2bw()`](https://ethz-ins.github.io/epiwraps/reference/bam2bw.md),
is not always appropriate. The main reasons are 1) that different
samples/experiments can have a different signal-to-noise ratio, with the
result that more sequencing is needed to obtain a similar coverage of
enriched region; 2) that there might be global differences in the amount
of the signal of interest (e.g. more or less binding, globally, in one
cell type vs another); and 3) that there might be differences in
technical biases, such as GC content. For these reasons, different
normalization methods are needed according to circumstances and what
assumptions seem reasonable. Here is an overview of the normalization
methods currently implemented in `epiwraps` via the
[`getNormFactors()`](https://ethz-ins.github.io/epiwraps/reference/getNormFactors.md)
function:

- The ‘background’ or ‘SES’ normalization method (they are synonyms
  here) assumes that the background noise should on average be the same
  across experiments ([Diaz et al., Stat Appl Gen et Mol Biol,
  2012](https://doi.org/10.1515/1544-6115.1750)), an assumption that
  works well in practice and is robust to global differences in the
  amount of signal when there are not large differences in
  signal-to-noise ratio.
- The ‘MAnorm’ approach ( [Shao et al., Genome Biology
  2012](https://doi.org/10.1186/gb-2012-13-3-r16)) assumes that regions
  that are commonly enriched (i.e.  common peaks) in two experiments
  should on average have the same signal in the two experiments.
- The ‘enriched’ approach assumes that enriched regions are on average
  similarly enriched across samples. Contrarily to ‘MAnorm’, these
  regions do not need to be in common across samples/experiments. This
  is not robust to global differences.
- The ‘top’ approach assumes that the maximum enrichment (after some
  trimming) in peaks is the same across samples/experiments.
- The ‘S3norm’ ([Xiang et al., NAR
  2020](https://doi.org/10.1093/nar/gkaa105)) and ‘2cLinear’ methods try
  to normalize both enrichment and background simultaneously. S3norm
  does this in a log-linear fashion (as in the publication), while
  ‘2cLinear’ does it on the original scale.

The normalization factors can be computed using
[`getNormFactors()`](https://ethz-ins.github.io/epiwraps/reference/getNormFactors.md)
:

``` r
suppressPackageStartupMessages(library(epiwraps))
# we fetch the path to the example bigwig file:
bwf <- system.file("extdata/example_atac.bw", package="epiwraps")
# we'll just double it to create a fake multi-sample dataset:
bwfiles <- c(atac1=bwf, atac2=bwf)
nf <- getNormFactors(bwfiles, method="background")
```

    ## Comparing coverage in random regions...

``` r
nf
```

    ## atac1 atac2 
    ##     1     1

In this case, since the files are identical, the factors are both 1.

Some normalization methods addtionally require peaks as input, e.g.:

``` r
peaks <- system.file("extdata/example_peaks.bed", package="epiwraps")
nf <- getNormFactors(bwfiles, peaks = peaks, method="MAnorm")
```

    ## Comparing coverage in peaks...

(Note that MAnorm would normally require to have a list of peaks for
each sample/experiment).

Once computed, the normalization factors can be applied to an
`EnrichmentSE` object:

``` r
sm <- signal2Matrix(bwfiles, peaks, extend=1000L)
```

    ## Reading /home/runner/work/_temp/Library/epiwraps/extdata/example_atac.bw
    ## Reading /home/runner/work/_temp/Library/epiwraps/extdata/example_atac.bw

``` r
sm <- renormalizeSignalMatrices(sm, scaleFactors=nf)
sm
```

    ## class: EnrichmentSE 
    ## 2 tracks across 150 regions
    ## assays(2): normalized input
    ## rownames(150): 1:195054101-195054250 1:133522798-133523047 ...
    ##   1:22224734-22224983 1:90375438-90375787
    ## rowData names(0):
    ## colnames(2): atac1 atac2
    ## colData names(0):
    ## metadata(0):

The object now has a new assay, called `normalized`, which has been put
in front and therefore will be used for most downstream usages unless
the uses specifies otherwise. Note that for any downstream function it
is however possible to specify which assay to use via the `assay`
argument.

## Obtaining normalization factors from the signal matrices themselves

It is also possible to normalize the signal matrices using factors
derived from the matrices themselves, using the
`renormalizeSignalMatrices` function. Note that this is provided as a
‘quick-and-dirty’ approach that does not have the robustness of proper
estimation methods. Specifically, beyond providing manual scaling
factors (e.g. computed using `getNormFactors`), the function includes
two methods :

- `method="border"` works on the assumption that the left/right borders
  of the matrices represent background signal which should be equal
  across samples. As such, it can be seen as an approximation of the
  aforementioned background normalization. However, it will work only
  if 1) the left/right borders of the matrices are sufficiently far from
  the signal (e.g. peaks) to be chiefly noise, and (as with the main
  background normalization method itself) 2) the signal-to-noise ratio
  is comparable across tracks/samples.
- `method="top"` instead works on the assumption that the highest signal
  (after some eventual trimming of the extremes) should be the same
  across tracks/samples.

To illustrate these, we will first introduce some difference between our
two tracks using arbitrary factors:

``` r
sm <- renormalizeSignalMatrices(sm, scaleFactors=c(1,4), toAssay="test")
plotEnrichedHeatmaps(sm, assay = "test")
```

![](normalization_files/figure-html/unnamed-chunk-5-1.png)

Then we can normalize:

``` r
sm <- renormalizeSignalMatrices(sm, method="top", fromAssay="test")
# again this adds an assay to the object, which will be automatically used when plotting:
plotEnrichedHeatmaps(sm)
```

    ## Using assay topNormalized

![](normalization_files/figure-html/unnamed-chunk-6-1.png)

And we’ve recovered comparable signal across the two tracks/samples.

  
  

## Session information

``` r
sessionInfo()
```

    ## R version 4.5.3 (2026-03-11)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] epiwraps_0.99.108           EnrichedHeatmap_1.40.1     
    ##  [3] ComplexHeatmap_2.26.1       SummarizedExperiment_1.40.0
    ##  [5] Biobase_2.70.0              GenomicRanges_1.62.1       
    ##  [7] Seqinfo_1.0.0               IRanges_2.44.0             
    ##  [9] S4Vectors_0.48.0            BiocGenerics_0.56.0        
    ## [11] generics_0.1.4              MatrixGenerics_1.22.0      
    ## [13] matrixStats_1.5.0           BiocStyle_2.38.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3          rstudioapi_0.18.0          
    ##   [3] jsonlite_2.0.0              shape_1.4.6.1              
    ##   [5] magrittr_2.0.4              magick_2.9.1               
    ##   [7] GenomicFeatures_1.62.0      farver_2.1.2               
    ##   [9] rmarkdown_2.31              GlobalOptions_0.1.3        
    ##  [11] fs_2.0.1                    BiocIO_1.20.0              
    ##  [13] ragg_1.5.2                  vctrs_0.7.2                
    ##  [15] memoise_2.0.1               Rsamtools_2.26.0           
    ##  [17] RCurl_1.98-1.18             base64enc_0.1-6            
    ##  [19] htmltools_0.5.9             S4Arrays_1.10.1            
    ##  [21] progress_1.2.3              curl_7.0.0                 
    ##  [23] SparseArray_1.10.10         Formula_1.2-5              
    ##  [25] sass_0.4.10                 bslib_0.10.0               
    ##  [27] htmlwidgets_1.6.4           desc_1.4.3                 
    ##  [29] Gviz_1.54.0                 httr2_1.2.2                
    ##  [31] cachem_1.1.0                GenomicAlignments_1.46.0   
    ##  [33] lifecycle_1.0.5             iterators_1.0.14           
    ##  [35] pkgconfig_2.0.3             Matrix_1.7-4               
    ##  [37] R6_2.6.1                    fastmap_1.2.0              
    ##  [39] clue_0.3-68                 digest_0.6.39              
    ##  [41] TFMPvalue_1.0.0             colorspace_2.1-2           
    ##  [43] AnnotationDbi_1.72.0        textshaping_1.0.5          
    ##  [45] Hmisc_5.2-5                 RSQLite_2.4.6              
    ##  [47] seqLogo_1.76.0              filelock_1.0.3             
    ##  [49] httr_1.4.8                  abind_1.4-8                
    ##  [51] compiler_4.5.3              bit64_4.6.0-1              
    ##  [53] doParallel_1.0.17           backports_1.5.0            
    ##  [55] htmlTable_2.4.3             S7_0.2.1                   
    ##  [57] BiocParallel_1.44.0         DBI_1.3.0                  
    ##  [59] biomaRt_2.66.2              rappdirs_0.3.4             
    ##  [61] DelayedArray_0.36.1         rjson_0.2.23               
    ##  [63] gtools_3.9.5                caTools_1.18.3             
    ##  [65] tools_4.5.3                 foreign_0.8-91             
    ##  [67] nnet_7.3-20                 glue_1.8.0                 
    ##  [69] restfulr_0.0.16             checkmate_2.3.4            
    ##  [71] cluster_2.1.8.2             TFBSTools_1.48.0           
    ##  [73] gtable_0.3.6                BSgenome_1.78.0            
    ##  [75] ensembldb_2.34.0            data.table_1.18.2.1        
    ##  [77] hms_1.1.4                   XVector_0.50.0             
    ##  [79] motifmatchr_1.32.0          foreach_1.5.2              
    ##  [81] pillar_1.11.1               stringr_1.6.0              
    ##  [83] limma_3.66.0                circlize_0.4.17            
    ##  [85] dplyr_1.2.0                 BiocFileCache_3.0.0        
    ##  [87] lattice_0.22-9              deldir_2.0-4               
    ##  [89] rtracklayer_1.70.1          bit_4.6.0                  
    ##  [91] biovizBase_1.58.0           DirichletMultinomial_1.52.0
    ##  [93] tidyselect_1.2.1            locfit_1.5-9.12            
    ##  [95] pbapply_1.7-4               Biostrings_2.78.0          
    ##  [97] knitr_1.51                  gridExtra_2.3              
    ##  [99] bookdown_0.46               ProtGenerics_1.42.0        
    ## [101] edgeR_4.8.2                 xfun_0.57                  
    ## [103] statmod_1.5.1               stringi_1.8.7              
    ## [105] UCSC.utils_1.6.1            lazyeval_0.2.2             
    ## [107] yaml_2.3.12                 evaluate_1.0.5             
    ## [109] codetools_0.2-20            cigarillo_1.0.0            
    ## [111] interp_1.1-6                GenomicFiles_1.46.0        
    ## [113] tibble_3.3.1                BiocManager_1.30.27        
    ## [115] cli_3.6.5                   rpart_4.1.24               
    ## [117] systemfonts_1.3.2           jquerylib_0.1.4            
    ## [119] dichromat_2.0-0.1           Rcpp_1.1.1                 
    ## [121] GenomeInfoDb_1.46.2         dbplyr_2.5.2               
    ## [123] png_0.1-9                   XML_3.99-0.23              
    ## [125] parallel_4.5.3              pkgdown_2.2.0              
    ## [127] ggplot2_4.0.2               blob_1.3.0                 
    ## [129] prettyunits_1.2.0           jpeg_0.1-11                
    ## [131] latticeExtra_0.6-31         AnnotationFilter_1.34.0    
    ## [133] bitops_1.0-9                pwalign_1.6.0              
    ## [135] viridisLite_0.4.3           VariantAnnotation_1.56.0   
    ## [137] scales_1.4.0                crayon_1.5.3               
    ## [139] GetoptLong_1.1.0            rlang_1.1.7                
    ## [141] cowplot_1.2.0               KEGGREST_1.50.0
