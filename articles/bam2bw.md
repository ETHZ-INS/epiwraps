# Generating bigwig tracks with bam2bw

Abstract

This vignette documents the use of the ‘bam2bw’ function of the epiwraps
package, which generates bigwig tracks from alignment files in an
efficient and flexible fashion.

## Introduction

The `bam2bw` function can be used to compute per-nucleotide or per-bin
coverage from alignments and save it to a bigwig file. In this process,
information about individual reads is lost, but the produced signals are
considerably more lightweight and amenable to visualization. The bigwig
format is readily queried from R or compatible with a variety of tools,
including genome browsers.

  

### Many ways of compiling coverage

To introduce the different variations on coverage, let’s assume you’ve
go the following single-end reads:

``` r
suppressPackageStartupMessages(library(epiwraps))
```

    ## Warning: replacing previous import 'IRanges::median' by 'stats::median' when
    ## loading 'epiwraps'

``` r
# we create some arbitrary genomic ranges
gr <- GRanges("chr1", IRanges(c(30,70,120), width=50), strand=c("+","+","-"),
              seqlengths=c(chr1=500))
plotSignalTracks(list(reads=gr), region="chr1:1:180", extend=0, genomeAxis=FALSE)
```

![](bam2bw_files/figure-html/unnamed-chunk-2-1.png)

For testing purposes, we’ll save this as a bam file and index it:

``` r
bam <- tempfile(fileext = ".bam") # temp file name
rtracklayer::export(gr, bam, format="bam")
Rsamtools::indexBam(bam)
```

    ##       /tmp/Rtmpqqcnw4/file3bd016d9c43c.bam 
    ## "/tmp/Rtmpqqcnw4/file3bd016d9c43c.bam.bai"

Using these example reads, we can illustrate different ways of computing
coverages.

First, we can save coverage at different resolutions, from full
resolution (each nucleotide is a single bin) to larger bin sizes, the
latter giving smaller filesizes:

``` r
# Full coverage with bin width of 1 nucleotide (i.e. full resolution)
cov_full_bw1 <- tempfile(fileext = ".bw") # temp file name
bam2bw(bam, cov_full_bw1, binWidth=1L, scaling=FALSE)
```

    ## `paired` not specified, assuming single-end reads. Set to paired='auto' to automatically detect.

    ## Reading in signal...

    ## Writing bigwig...

``` r
# Full coverage with larger bins
cov_full_bw25 <- tempfile(fileext = ".bw")
bam2bw(bam, cov_full_bw25, binWidth=25L, scaling=FALSE)
```

    ## `paired` not specified, assuming single-end reads. Set to paired='auto' to automatically detect.

    ## Reading in signal...

    ## Writing bigwig...

``` r
plotSignalTracks(list(reads=gr, "binWidth=1"=cov_full_bw1, "binWidth=25"=cov_full_bw25),
                 region="chr1:1:180", extend=0)
```

![](bam2bw_files/figure-html/unnamed-chunk-4-1.png)

Both tracks compile the number of reads that overlap each position, but
in the bottom track the signal is by chunks of 25 nucleotides. By
default, the maximum signal inside a bin is used, however it is possible
to change this:

``` r
# Using mean per bin:
cov_full_bw25mean <- tempfile(fileext = ".bw")
bam2bw(bam, cov_full_bw25mean, binWidth=25L, binSummarization = "mean", scaling=FALSE)
```

    ## `paired` not specified, assuming single-end reads. Set to paired='auto' to automatically detect.

    ## Reading in signal...

    ## Writing bigwig...

``` r
plotSignalTracks(list(reads=gr, "binWidth=1"=cov_full_bw1, 
                      "binWidth=25\n(max)"=cov_full_bw25,
                      "binWidth=25\n(mean)"=cov_full_bw25mean),
                 region="chr1:1:180", extend=0)
```

![](bam2bw_files/figure-html/unnamed-chunk-5-1.png)

In most cases, single-end reads are just the beginning of the DNA
fragments obtained, and we know the average size of the fragments from
library prep QC (if not, see the `estimateFragSize` function, or the
simpler `estimate.mean.fraglen` function of the
*[chipseq](https://bioconductor.org/packages/3.22/chipseq)* package). It
is therefore common to extend reads to this size when computing
coverage, so as to obtain the number of fragments (rather than reads)
coverage each position. This can be done as follows:

``` r
# Here the reads are 50bp, and we want to extend them to 100bp, hence _by_ 50:
cov_full_ext <- tempfile(fileext = ".bw")
bam2bw(bam, cov_full_ext, binWidth=1L, extend=50L, scaling=FALSE)
```

    ## `paired` not specified, assuming single-end reads. Set to paired='auto' to automatically detect.

    ## Reading in signal...

    ## Writing bigwig...

``` r
plotSignalTracks(list(reads=gr, "no extension"=cov_full_bw1, 
                      "read extension"=cov_full_ext),
                 region="chr1:1:190", extend=0)
```

![](bam2bw_files/figure-html/unnamed-chunk-6-1.png)

Taking into account read extension shows better the peak at the center,
which is indicative of the fact that this is the region that was the
most enriched in the captured fragments.

Instead of computing coverage, we could compute the number of reads
starting, ending, or being centered at each position:

``` r
# Here the reads are 50bp, and we want to extend them to 100bp, hence _by_ 50:
cov_start <- tempfile(fileext = ".bw")
bam2bw(bam, cov_start, binWidth=1L, extend=50L, scaling=FALSE, type="start")
```

    ## `paired` not specified, assuming single-end reads. Set to paired='auto' to automatically detect.

    ## Reading in signal...

    ## Writing bigwig...

``` r
cov_center <- tempfile(fileext = ".bw")
bam2bw(bam, cov_center, binWidth=1L, extend=50L, scaling=FALSE, type="center")
```

    ## `paired` not specified, assuming single-end reads. Set to paired='auto' to automatically detect.

    ## Reading in signal...

    ## Writing bigwig...

``` r
plotSignalTracks(list(reads=gr, "type=full"=cov_full_bw1, 
                      "type=start"=cov_start, "type=center"=cov_center),
                 region="chr1:1:190", extend=0)
```

![](bam2bw_files/figure-html/unnamed-chunk-7-1.png)

Note that when extending reads, as in this case, the position
(e.g. “center”) are relative to the extended read (i.e. extension is
applied first).

  

### Example heatmaps created using different bigwig generation procedures

The following figure, created using `epiwraps` (see [the vignette on
generating such
plots](https://ethz-ins.github.io/epiwraps/articles/multiRegionPlot.md)),
represent chromatin accessibility (ATAC-seq) signals around bound CTCF
motifs in T-cells. The different signals are based on different bigwig
files derived from the same bam file using the functions described
above.

![](example_bigwigs_hm.png)

- The first heatmap (‘full coverage’) was generated with default
  parameter, and is the fragment coverage (i.e. how many fragments
  overlap any given location).
- The second heatmap shows the fragment of sizes compatible with
  mono-nucleosomes, resizing fragments from their centers. The exact
  arguments used were
  `minFragLength=147, maxFragLength=230, type="center", extend=25L`.
  Using this we can see nucleosomes well-positioned at some distance
  from CTCF binding sites, but a relative depletion of nucleosomes at
  the bound site itself.
- The third shows the coverage of nucleosome-free fragments, in this
  case it was used with `maxFragLength=120`. These are indicative of TF
  binding, and indeed we only see and enrichment in the center.
- The fourth shows where the transposase inserted itself, and was
  generated with `trim=4L, binWidth=1L, maxFragLength=120, type="ends"`.
  Using this we can see a nice footprint protected from the transposase
  by CTCF binding.

For more information about these parameters, see
[`?bam2bw`](https://ethz-ins.github.io/epiwraps/reference/bam2bw.md).

ATAC fragment sizes convey rich information, but they are not perfect.
For example, regions of DNA bound by multiple TFs can sometimes be
captured as longer fragments which we would mistakenly classify as
mono-nucleosome containing.

  

## Working with fragment files as an input

If you use fragment files (preferably tabix-indexed) rather than bam
files as input, you can still perform most of the above tasks. See the
[`?frag2bw`](https://ethz-ins.github.io/epiwraps/reference/frag2bw.md)
function for more information.

  
  

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
