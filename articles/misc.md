# Miscellaneous epiwraps functions

Abstract

This vignette introduces some more isolated epiwraps functions.

## Quality control

### Coverage statistics

Coverage statistics give an overview of how the reads are distributed
across the genome (or more precisely, across a large number of random
regions). The `getCovStats` will compute such statistics from bam or
bigwig files (from bigwig files will be considerably faster, but if the
files are normalized the coverage/density will be relative).

Because our example data spans only part of a chromosome, we’ll exclude
completely empty regions using the `exclude` parameter, which would
normally be used to exclude regions likely to be technical artefacts
(e.g. blacklisted regions).

``` r
suppressPackageStartupMessages(library(epiwraps))
# get the path to an example bigwig file:
bwf <- system.file("extdata/example_atac.bw", package="epiwraps")
cs <- getCovStats(bwf, exclude=GRanges("1", IRanges(1, 4300000)))
plotCovStats(cs)
```

![](misc_files/figure-html/covstats-1.png)

Panel A shows the proportion of sampled regions which are above a
certain read density (relative because this is a normalized bigwig file,
would be coverage otherwise). This shows us, for example, that as
expected only a minority of regions have any reads at all (indicating
that the reads are not randomly distributed). Panel B is what is
sometimes called a fingerprint plot. It similarly shows us that the
reads are concentrated in very few regions, since the vast majority of
regions have only a very low fraction of the coverage of a few
high-density regions. Randomly distributed reads would go along the
diagonal, but one normally has a curve somewhere between this line and
the lower-right corner – the farther away from the diagonal, to more
strongly enriched the data is.

This can be done for multiple files simultaneously. If we have several
files, we can also use the coverage in the random windows to computer
their similarity (see
[`?plotCorFromCovStats`](https://ethz-ins.github.io/epiwraps/reference/plotCorFromCovStats.md)).

### Fragment length distributions

Given one or more paired-end bam files, we can extract and plot the
fragment length distribution using:

``` r
fragSizesDist(bam, what=100)
```

### TSS enrichment

The TSS enrichment can also be calculated and plotted using the
`TSSenrichment` function.

## Peak calling

A very experimental peak calling function can be used, either against an
input control or against local or global backgrounds:

``` r
p <- callPeaks(bam, fragLength=50)
```

Note that this function is still under heavy development, and its usage
at the moment is discouraged!

## Region overlapping

The `GenomicRanges` package offers fast and powerful functions for
overlapping genomic regions. `epiwraps` includes wrappers around those
for common tasks, such as calculating and visualizing overlaps across
multiple sets of regions (see `?regionUpset`,
[`?regionOverlaps`](https://ethz-ins.github.io/epiwraps/reference/regionOverlaps.md),
and
[`?regionCAT`](https://ethz-ins.github.io/epiwraps/reference/regionCAT.md)).

  
  

## Session information

``` r
sessionInfo()
```

    ## R version 4.5.3 (2026-03-11)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
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
    ##  [1] epiwraps_0.99.106           EnrichedHeatmap_1.40.1     
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
    ##   [5] magrittr_2.0.4              GenomicFeatures_1.62.0     
    ##   [7] farver_2.1.2                rmarkdown_2.30             
    ##   [9] GlobalOptions_0.1.3         fs_2.0.0                   
    ##  [11] BiocIO_1.20.0               ragg_1.5.1                 
    ##  [13] vctrs_0.7.2                 memoise_2.0.1              
    ##  [15] Rsamtools_2.26.0            RCurl_1.98-1.18            
    ##  [17] base64enc_0.1-6             htmltools_0.5.9            
    ##  [19] S4Arrays_1.10.1             progress_1.2.3             
    ##  [21] curl_7.0.0                  SparseArray_1.10.9         
    ##  [23] Formula_1.2-5               sass_0.4.10                
    ##  [25] bslib_0.10.0                htmlwidgets_1.6.4          
    ##  [27] desc_1.4.3                  Gviz_1.54.0                
    ##  [29] httr2_1.2.2                 cachem_1.1.0               
    ##  [31] GenomicAlignments_1.46.0    lifecycle_1.0.5            
    ##  [33] iterators_1.0.14            pkgconfig_2.0.3            
    ##  [35] Matrix_1.7-4                R6_2.6.1                   
    ##  [37] fastmap_1.2.0               clue_0.3-67                
    ##  [39] digest_0.6.39               TFMPvalue_1.0.0            
    ##  [41] colorspace_2.1-2            AnnotationDbi_1.72.0       
    ##  [43] textshaping_1.0.5           Hmisc_5.2-5                
    ##  [45] RSQLite_2.4.6               labeling_0.4.3             
    ##  [47] seqLogo_1.76.0              filelock_1.0.3             
    ##  [49] httr_1.4.8                  abind_1.4-8                
    ##  [51] compiler_4.5.3              withr_3.0.2                
    ##  [53] bit64_4.6.0-1               doParallel_1.0.17          
    ##  [55] backports_1.5.0             htmlTable_2.4.3            
    ##  [57] S7_0.2.1                    BiocParallel_1.44.0        
    ##  [59] DBI_1.3.0                   biomaRt_2.66.2             
    ##  [61] rappdirs_0.3.4              DelayedArray_0.36.0        
    ##  [63] rjson_0.2.23                gtools_3.9.5               
    ##  [65] caTools_1.18.3              tools_4.5.3                
    ##  [67] foreign_0.8-91              nnet_7.3-20                
    ##  [69] glue_1.8.0                  restfulr_0.0.16            
    ##  [71] checkmate_2.3.4             cluster_2.1.8.2            
    ##  [73] TFBSTools_1.48.0            gtable_0.3.6               
    ##  [75] BSgenome_1.78.0             ensembldb_2.34.0           
    ##  [77] data.table_1.18.2.1         hms_1.1.4                  
    ##  [79] XVector_0.50.0              motifmatchr_1.32.0         
    ##  [81] foreach_1.5.2               pillar_1.11.1              
    ##  [83] stringr_1.6.0               circlize_0.4.17            
    ##  [85] dplyr_1.2.0                 BiocFileCache_3.0.0        
    ##  [87] lattice_0.22-9              deldir_2.0-4               
    ##  [89] rtracklayer_1.70.1          bit_4.6.0                  
    ##  [91] biovizBase_1.58.0           DirichletMultinomial_1.52.0
    ##  [93] tidyselect_1.2.1            locfit_1.5-9.12            
    ##  [95] pbapply_1.7-4               Biostrings_2.78.0          
    ##  [97] knitr_1.51                  gridExtra_2.3              
    ##  [99] bookdown_0.46               ProtGenerics_1.42.0        
    ## [101] xfun_0.57                   stringi_1.8.7              
    ## [103] UCSC.utils_1.6.1            lazyeval_0.2.2             
    ## [105] yaml_2.3.12                 evaluate_1.0.5             
    ## [107] codetools_0.2-20            cigarillo_1.0.0            
    ## [109] interp_1.1-6                GenomicFiles_1.46.0        
    ## [111] tibble_3.3.1                BiocManager_1.30.27        
    ## [113] cli_3.6.5                   rpart_4.1.24               
    ## [115] systemfonts_1.3.2           jquerylib_0.1.4            
    ## [117] dichromat_2.0-0.1           Rcpp_1.1.1                 
    ## [119] GenomeInfoDb_1.46.2         dbplyr_2.5.2               
    ## [121] png_0.1-9                   XML_3.99-0.23              
    ## [123] parallel_4.5.3              pkgdown_2.2.0              
    ## [125] ggplot2_4.0.2               blob_1.3.0                 
    ## [127] prettyunits_1.2.0           jpeg_0.1-11                
    ## [129] latticeExtra_0.6-31         AnnotationFilter_1.34.0    
    ## [131] bitops_1.0-9                pwalign_1.6.0              
    ## [133] viridisLite_0.4.3           VariantAnnotation_1.56.0   
    ## [135] scales_1.4.0                crayon_1.5.3               
    ## [137] GetoptLong_1.1.0            rlang_1.1.7                
    ## [139] cowplot_1.2.0               KEGGREST_1.50.0
