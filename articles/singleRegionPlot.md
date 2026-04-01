# Visualizing signals in a single region

Abstract

This vignette documents the use of the ‘plotSignalTracks’ to generate
genome-browser-like plots of signals and annotations along genomic
coordinates in a single given region. It is chiefly a wrapper around the
‘Gviz’ package.

## Plotting signals in a region

The `plotSignalTracks` function is a wrapper around the
*[Gviz](https://bioconductor.org/packages/3.22/Gviz)* package, which
plots one or more signals along genomic coordinates (in a genome-browser
like fashion). The function lacks the full flexibility of the
*[Gviz](https://bioconductor.org/packages/3.22/Gviz)* package, but
presents a considerable simpler interface, with automatic default
parameters, etc. It has two essential arguments: a (named) list of files
whose signal to display (can be a mixture of bigwig, bam, or bed-like
files), and the region in which to display the signals (can be given as
a GRanges or as a string). The function then automatically determines
the relevant track type and setting from the file types.

``` r
suppressPackageStartupMessages(library(epiwraps))
```

    ## Warning: replacing previous import 'IRanges::median' by 'stats::median' when
    ## loading 'epiwraps'

``` r
# get the path to an example bigwig file:
bwf1 <- system.file("extdata/example_rna.bw", package="epiwraps")
plotSignalTracks(list(RNA=bwf1), region="8:22165140-22212326", genomeAxis=TRUE)
```

![](singleRegionPlot_files/figure-html/signal1-1.png)

``` r
# we could plot multiple tracks as follows:
plotSignalTracks(list(track1=bwf1, track2=bwf1), region="8:22165140-22212326")
```

![](singleRegionPlot_files/figure-html/signal1-2.png)

`GRanges` objects can also be plotted as annotation tracks alongside
other data:

``` r
myregions <- GRanges("8", IRanges(c(22166000,22202300), width=3000))
plotSignalTracks(list(RNA=bwf1, regions=myregions), region="8:22165140-22212326")
```

![](singleRegionPlot_files/figure-html/signal2-1.png)

Colors, track display types, and such parameters can either be set for
all tracks or for each individual track, for example:

``` r
myregions <- GRanges("8", IRanges(c(22166000,22202300), width=3000))
plotSignalTracks(list(RNA=bwf1, regions=myregions), colors=c("red", "black"),
                 region="8:22165140-22212326")
```

![](singleRegionPlot_files/figure-html/signal3-1.png)

For bam files, we can also plot individual reads:

``` r
# we fetch an example bam file:
bam <- system.file("extdata", "ex1.bam", package="Rsamtools")
plotSignalTracks(c("my bam file"=bam), "seq1:1-1500", type="alignments")
```

    ## Warning in call_new_fun_in_cigarillo("sequenceLayer", "project_sequences", : sequenceLayer() is formally deprecated in GenomicAlignments >= 1.45.5 and
    ##   replaced with the project_sequences() function from the new cigarillo package

![](singleRegionPlot_files/figure-html/signalAlignments-1.png)

### Merging signal from different tracks

In addition to being displayed one below the other, data tracks can be
combined in different ways. To do this, the tracks can simply be given
in a nested fashion:

``` r
plotSignalTracks(list(track1=bwf1, combined=c(bwf1,bwf1)),
                 region="8:22165140-22212326")
```

In this example we are always using the same track, but the first
element (‘track1’) plots the track alone, while the second (‘combined’)
merges the two given tracks. By default, the mean is shown, but this can
be controlled through the `aggregation` argument. In addition to usual
operations, the tracks can be overlayed on top of one another
(`aggregation='overlay'`), or shown as a heatmap
(`aggregation='heatmap'`), or a mean track with heatmap below
(`aggregation="heatmap+mean"`).

### Using an EnsDb object

If an `EnsDb` object is available (see the
*[ensembldb](https://bioconductor.org/packages/3.22/ensembldb)* package
for a description of the class and its methods, and the
*[AnnotationHub](https://bioconductor.org/packages/3.22/AnnotationHub)*
package for a convenient way of fetching such annotation objects), two
additional options are available: first, instead of specifying the
region as coordinates, one can specify a gene or transcript name, and
the corresponding region will be fetched. In addition, the genes or
transcripts can be displayed. For example:

``` r
# we fetch the GRCh38 Ensembl 103 annotation (this is not run in the vignette,
# as it takes some time to download the annotation the first time is used):
library(AnnotationHub)
ah <- AnnotationHub()
ensdb <- ah[["AH89426"]]
# we plot our previous RNA bigwig file, around the BMP1 locus:
plotSignalTracks(c(coverage=bwf1), region="BMP1", ensdb=ensdb, 
                 transcripts="full")
```

![](TracksWithTranscripts.png)

Now we can see that the coverage is nicely restricted to exons, and that
some transcripts/exons are not expressed as highly as others. The
transcripts could also have been collapsed into a gene model using
`transcripts="collapsed"` (the default).

To display only the gene track, the first argument can simply be
omitted.

### Further track customization

In addition to the `colors` and `type` argument (and a number of
others), which can customize the appearance of tracks, any additional
parameters supported by the respective
*[Gviz](https://bioconductor.org/packages/3.22/Gviz)* function can be
passed through the `genes.params` (for Gviz’s `GeneRegionTrack`),
`align.params` (for Gviz’s `AlignmentsTrack`, when plotting individual
reads), or `tracks.params` (for any other Gviz `DataTrack`).

For example, if you wish to manually set the same y-axis range for all
data tracks, this can be done with:

``` r
plotSignalTracks(list(track1=bwf1, track2=bwf1), region="8:22165140-22212326",
                 tracks.params=list(ylim=c(0,200)))
```

![](singleRegionPlot_files/figure-html/unnamed-chunk-2-1.png)

Also, in addition to passing filepaths or `GRanges`, any Gviz track(s)
can be passed (i.e. objects inheriting the `GdObject` class) can be
passed, enabling full track customization when needed.

  
  

## Session info

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
    ##   [5] magrittr_2.0.4              GenomicFeatures_1.62.0     
    ##   [7] farver_2.1.2                rmarkdown_2.31             
    ##   [9] GlobalOptions_0.1.3         fs_2.0.1                   
    ##  [11] BiocIO_1.20.0               ragg_1.5.2                 
    ##  [13] vctrs_0.7.2                 memoise_2.0.1              
    ##  [15] Rsamtools_2.26.0            RCurl_1.98-1.18            
    ##  [17] base64enc_0.1-6             htmltools_0.5.9            
    ##  [19] S4Arrays_1.10.1             progress_1.2.3             
    ##  [21] curl_7.0.0                  SparseArray_1.10.10        
    ##  [23] Formula_1.2-5               sass_0.4.10                
    ##  [25] bslib_0.10.0                htmlwidgets_1.6.4          
    ##  [27] desc_1.4.3                  Gviz_1.54.0                
    ##  [29] httr2_1.2.2                 cachem_1.1.0               
    ##  [31] GenomicAlignments_1.46.0    lifecycle_1.0.5            
    ##  [33] iterators_1.0.14            pkgconfig_2.0.3            
    ##  [35] Matrix_1.7-4                R6_2.6.1                   
    ##  [37] fastmap_1.2.0               clue_0.3-68                
    ##  [39] digest_0.6.39               TFMPvalue_1.0.0            
    ##  [41] colorspace_2.1-2            AnnotationDbi_1.72.0       
    ##  [43] textshaping_1.0.5           Hmisc_5.2-5                
    ##  [45] RSQLite_2.4.6               seqLogo_1.76.0             
    ##  [47] filelock_1.0.3              httr_1.4.8                 
    ##  [49] abind_1.4-8                 compiler_4.5.3             
    ##  [51] bit64_4.6.0-1               doParallel_1.0.17          
    ##  [53] backports_1.5.0             htmlTable_2.4.3            
    ##  [55] S7_0.2.1                    BiocParallel_1.44.0        
    ##  [57] DBI_1.3.0                   biomaRt_2.66.2             
    ##  [59] rappdirs_0.3.4              DelayedArray_0.36.1        
    ##  [61] rjson_0.2.23                gtools_3.9.5               
    ##  [63] caTools_1.18.3              tools_4.5.3                
    ##  [65] foreign_0.8-91              nnet_7.3-20                
    ##  [67] glue_1.8.0                  restfulr_0.0.16            
    ##  [69] checkmate_2.3.4             cluster_2.1.8.2            
    ##  [71] TFBSTools_1.48.0            gtable_0.3.6               
    ##  [73] BSgenome_1.78.0             ensembldb_2.34.0           
    ##  [75] data.table_1.18.2.1         hms_1.1.4                  
    ##  [77] XVector_0.50.0              motifmatchr_1.32.0         
    ##  [79] foreach_1.5.2               pillar_1.11.1              
    ##  [81] stringr_1.6.0               circlize_0.4.17            
    ##  [83] dplyr_1.2.0                 BiocFileCache_3.0.0        
    ##  [85] lattice_0.22-9              deldir_2.0-4               
    ##  [87] rtracklayer_1.70.1          bit_4.6.0                  
    ##  [89] biovizBase_1.58.0           DirichletMultinomial_1.52.0
    ##  [91] tidyselect_1.2.1            locfit_1.5-9.12            
    ##  [93] pbapply_1.7-4               Biostrings_2.78.0          
    ##  [95] knitr_1.51                  gridExtra_2.3              
    ##  [97] bookdown_0.46               ProtGenerics_1.42.0        
    ##  [99] xfun_0.57                   stringi_1.8.7              
    ## [101] UCSC.utils_1.6.1            lazyeval_0.2.2             
    ## [103] yaml_2.3.12                 evaluate_1.0.5             
    ## [105] codetools_0.2-20            cigarillo_1.0.0            
    ## [107] interp_1.1-6                GenomicFiles_1.46.0        
    ## [109] tibble_3.3.1                BiocManager_1.30.27        
    ## [111] cli_3.6.5                   rpart_4.1.24               
    ## [113] systemfonts_1.3.2           jquerylib_0.1.4            
    ## [115] dichromat_2.0-0.1           Rcpp_1.1.1                 
    ## [117] GenomeInfoDb_1.46.2         dbplyr_2.5.2               
    ## [119] png_0.1-9                   XML_3.99-0.23              
    ## [121] parallel_4.5.3              pkgdown_2.2.0              
    ## [123] ggplot2_4.0.2               blob_1.3.0                 
    ## [125] prettyunits_1.2.0           jpeg_0.1-11                
    ## [127] latticeExtra_0.6-31         AnnotationFilter_1.34.0    
    ## [129] bitops_1.0-9                pwalign_1.6.0              
    ## [131] viridisLite_0.4.3           VariantAnnotation_1.56.0   
    ## [133] scales_1.4.0                crayon_1.5.3               
    ## [135] GetoptLong_1.1.0            rlang_1.1.7                
    ## [137] cowplot_1.2.0               KEGGREST_1.50.0
