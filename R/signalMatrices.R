#' signal2Matrix
#' 
#' Reads the signals around (the centers of) a set of regions.
#'
#' @param filepaths A named vector of filepaths (preferentially to bigwig files)
#' @param regions A `GRanges` of the regions/positions around which to plot, or
#' the path to a bed file of such regions.
#' @param extend Number of basepair to extend on either side of the regions
#' @param w Bin size
#' @param cuts Whether to count cuts (e.g. beginning/end of fragments) rather
#' than coverage (ignored unless the input are bam files)
#' @param RPM Whether to perform RPM normalization (for bam input)
#' @param nthreads The number of threads to use to read and prepare the data 
#' (default 1). Alternatively, a \code{\link[BiocParallel]{BiocParallelParam}}
#' object.
#' @param flgs Flags for bam reading
#' @param ... Passed to \code{\link[EnrichedHeatmap]{normalizeToMatrix}}.
#'
#' @return A list of `normalizeToMatrix` objects
#' @export
#' 
#' @import GenomicRanges
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam
#' @importFrom Rsamtools scanBamFlag ScanBamParam countBam
#' @importFrom genomation ScoreMatrix
#' @import EnrichedHeatmap
#' @importFrom rtracklayer import BigWigSelection
#' @importFrom GenomicAlignments readGAlignmentPairs
#' 
#' @examples 
#' # we fetch the path to the example bigwig file:
#' (bw <- system.file("extdata/example_atac.bw", package="epiwraps"))
#' # we load example regions:
#' regions <- rtracklayer::import(system.file("extdata/example_peaks.bed", 
#'                                            package="epiwraps"))
#' length(regions)
#' # we obtain the matrix of the signal around the regions:
#' m <- signal2Matrix(bw, regions)
#' dim(m[[1]])
#' # we can plot it with:
#' plotEnrichedHeatmaps(m)
#' # we could also take a broader range around the center of the regions, and 
#' # use bigger bins:
#' m <- signal2Matrix(bw, regions, extend=2000, w=20)
#' # the matrix has the same size, but shows broader regions:
#' dim(m[[1]])
#' plotEnrichedHeatmaps(m)
signal2Matrix <- function(filepaths, regions, extend=1000, w=10, cuts=FALSE, 
                          ..., RPM=TRUE, nthreads=1L,
                          flgs=scanBamFlag(isDuplicate=FALSE,
                                           isSecondaryAlignment=FALSE)){
  if(cuts && !all(grepl("\\.bam$",filepaths,ignore.case=TRUE)))
    stop("`cuts` can only be used with BAM files.")
  if(is.null(names(filepaths)))
    names(filepaths) <- .cleanFileNames(filepaths)
  
  if(is.character(regions)){
    stopifnot(is.character(regions) && length(regions)==1)
    regions <- import(regions)
  }
  stopifnot(is(regions, "GRanges"))
  regions2 <- resize(regions, fix="center", width=extend*2)
  
  if(is.numeric(nthreads)){
    if(nthreads<2L){
      nthreads <- BiocParallel::SerialParam()
    }else{
      nthreads <- BiocParallel::MulticoreParam(as.integer(nthreads))
    }
  }
  
  bplapply(setNames(names(filepaths),names(filepaths)), BPPARAM=nthreads, 
           FUN=function(filename){
    filepath <- filepaths[[filename]]
    message("Reading ", filepath)
    if(grepl("\\.bam$",filepath,ignore.case=TRUE)){
      ####### BAM INPUT
      
      libsize <- NULL
      if(is.null(RPM) || RPM)
        libsize <- Rsamtools::countBam(filepath, 
                                       param=ScanBamParam(flag=flgs))$records
      if(cuts){
        params <- Rsamtools::ScanBamParam(which=regions2, flag=flgs)
        bam <- GenomicAlignments::readGAlignmentPairs(filepath, param=params)
        bam <- coverage(.align2cuts(bam))
        mat <- genomation::ScoreMatrix(bam, windows=regions2, 
                                       strand.aware=FALSE, library.size=libsize)
        rm(bam)
        if(is.null(RPM) || RPM) mat <- mat*1000000/libsize
      }else{
        mat <- genomation::ScoreMatrix(filepath, regions2, bam.paired.end=TRUE, 
                                       unique=TRUE, library.size=libsize)
        if(is.null(RPM) || RPM) mat <- mat*1000000/libsize
      }
      if(w!=1){
        nwin <- round((2*extend)/w)
        if(nwin %% 2 != 0) nwin <- nwin+1
        mat <- sapply(split(seq_len(ncol(mat)),cut(seq_len(ncol(mat)),nwin)), 
                      FUN=function(x) rowMeans(mat[,x]) )
        extend <- floor(ncol(mat)/2)
      }
      mat <- EnrichedHeatmap::as.normalizedMatrix( 
        unclass(mat), extend=extend, signal_name=filename, k_target=0, 
        k_upstream=extend, k_downstream=extend+(w==1) )
      
      ####### END BAM INPUT
      
    }else if(grepl("\\.bw$|\\.bigwig$", filepath)){
      
      ####### BIGWIG INPUT
      gr <- rtracklayer::import(filepath, format="BigWig", 
                                selection=BigWigSelection(regions2))
      mat <- EnrichedHeatmap::normalizeToMatrix( signal=gr, ..., w=w,
                                target=resize(regions, fix="center", width=1),
                                mean_mode="w0", target_ratio = 0,
                                value_column="score", extend=extend )
      ####### END BIGWIG INPUT
      
    }else{
      stop("Unknown file format")
    }
    
    mat
  })
}


#' plotEnrichedHeatmaps
#' 
#' Plots enrichment heatmaps from the output of `signal2Matrix`
#'
#' @param ml A named list of matrices as produced by \code{\link{signal2Matrix}}
#' @param trim The quantile above which to trim values for the colorscale
#' @param colors The heatmap colors to use
#' @param scale_title The title of the scale
#' @param title_size The size of the heatmap titles
#' @param row_order Option order of the rows
#' @param cluster_rows Whether to cluster rows
#' @param ... Passed to EnrichedHeatmap::EnrichedHeatmap
#'
#' @return A HeatmapList object
#' @import EnrichedHeatmap
#' @importFrom viridisLite inferno
#' @importFrom circlize colorRamp2
#' @examples 
#' # we first fetch the path to the example bigwig file:
#' bw <- system.file("extdata/example_atac.bw", package="epiwraps")
#' # Since we only have one, we'll use the same and pretend they're 2 samples:
#' bw <- c(sample1=bw, sample2=bw)
#' # we next load regions of interest (either GRanges or path to a bed file):
#' regions <- system.file("extdata/example_peaks.bed", package="epiwraps")
#' # we obtain the matrix of the signal around the regions:
#' m <- signal2Matrix(bw, regions)
#' plotEnrichedHeatmaps(m)
#' # we could also just plot one with:
#' # plotEnrichedHeatmaps(m[1])
#' # or change the aesthetics, e.g.:
#' plotEnrichedHeatmaps(m, trim=0.98, scale_title="RPKM", 
#'                      colors=c("white","darkred"))
#' # any argument accepted by `EnrichedHeatmap` (and hence by 
#' # `ComplexHeatmap::Heatmap`) can be used, e.g.: 
#' plotEnrichedHeatmaps(m, row_title="My regions of interest")
plotEnrichedHeatmaps <- function(ml, trim=0.998, colors=inferno(100), 
                                 scale_title="density", title_size=11, 
                                 row_order=NULL, cluster_rows=FALSE, ...){
  if(!is.list(ml)) ml <- list(signal=ml)
  hl <- NULL
  ylim <- c(0,max(unlist(lapply(ml,FUN=function(x) max(colMeans(x))))))
  common_min <- min(unlist(lapply(ml,min)))
  common_max <- max(unlist(lapply(ml,prob=trim,FUN=quantile)))
  breaks <- seq(from=common_min, to=common_max, length.out=length(colors))
  col_fun <- circlize::colorRamp2(breaks, colors)
  if(isTRUE(cluster_rows)){
    cluster_rows <- hclust(dist(do.call(cbind, ml)))
  }else if(is.null(row_order)){
    row_order <- order(-rowMeans(do.call(cbind, lapply(ml, enriched_score))))
  }
  for(m in names(ml)){
    isLast <- m==rev(names(ml))[1]
    HA <- HeatmapAnnotation(enriched=anno_enriched(ylim=ylim, show_error=TRUE,
                                                   axis=isLast))
    hl <- hl + EnrichedHeatmap(ml[[m]], column_title=m, col=col_fun, ...,
                    column_title_gp=gpar(fontsize=title_size),
                    cluster_rows = cluster_rows, row_order=row_order,
                    top_annotation=HA, show_heatmap_legend=isLast,
                    name=ifelse(isLast,scale_title,m) )
  }
  hl
}


#' meltSignals
#'
#' Aggregates and melts a list of signal matrices, for plotting (with ggplot).
#'
#' @param ml A named list of matrices as produced by \code{\link{signal2Matrix}}
#' @param fun The aggregation to perform. Either "mean", "sum", "median", or
#' a custom function to be applied on columns.
#'
#' @return A data.frame.
#' @export
#' @importFrom matrixStats colMedians
meltSignals <- function(ml, fun=c("mean","sum","median")){
  stopifnot(is.list(ml))
  if(!is.function(fun)) fun <- match.arg(fun)
  a <- attributes(ml[[1]])
  ai <- c("upsteam_index", "downsteam_index", "target_index", "extend")
  if(!all(unlist(lapply(ml, FUN=function(x){
    identical(attributes(x)[ai],a[ai]) && nrow(x)==nrow(ml[[1]])
  })))) stop("The matrices do not have an homogeneous design!")
  x <- c( -1*(a$extend[1]/length(a$upstream_index))*
            rev(seq_along(a$upstream_index)),
          (a$extend[2]/length(a$downstream_index))*
            seq_along(a$downstream_index) )
  y <- lapply(ml, FUN=function(x){
    if(is.function(fun)) return(apply(x,2,fun))
    switch(fun,
           mean=colMeans(x),
           sum=colSums(x),
           median=matrixStats::colMedians(x))
  })
  d <- data.frame( position=rep(x, length(ml)), value=unlist(y),
                   sample=rep(factor(names(ml), names(ml)), each=length(x)) )
  d
}

