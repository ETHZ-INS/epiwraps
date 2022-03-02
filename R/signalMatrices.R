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
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object, or 
#' the number of threads to use to read and prepare the data.
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
                          ..., RPM=TRUE, BPPARAM=1L,
                          flgs=scanBamFlag(isDuplicate=FALSE,
                                           isSecondaryAlignment=FALSE)){
  if(!all(unlist(lapply(filepaths, file.exists))))
    stop("Some of the files given do not exist, check the paths provided.")
  if(cuts && !all(grepl("\\.bam$",filepaths,ignore.case=TRUE)))
    stop("`cuts` can only be used with BAM files.")
  if(is.null(names(filepaths)))
    names(filepaths) <- .cleanFileNames(filepaths)
  
  if(is.character(regions)){
    stopifnot(is.character(regions) && length(regions)==1)
    if(!file.exists(regions)) stop("The `regions` file appears not to exist.")
    regions <- import(regions)
  }
  stopifnot(is(regions, "GRanges"))
  if(is.null(names(regions)))
    names(regions) <- paste0("region", seq_along(regions))
  regions2 <- resize(regions, fix="center", width=extend*2)
  
  bplapply(setNames(names(filepaths),names(filepaths)), BPPARAM=.getBP(BPPARAM), 
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
    tryCatch({
      names(mat) <- names(regions)
      mat
    }, error=function(e){
      warning(e)
      mat
    })
  })
}


#' plotEnrichedHeatmaps
#' 
#' Plots enrichment heatmaps from the output of `signal2Matrix`. This is a 
#' convenience wrapper around \code{\link[EnrichedHeatmap]{EnrichedHeatmap}}.
#'
#' @param ml A named matrix list as produced by \code{\link{signal2Matrix}}.
#' @param trim The quantile above which to trim values for the colorscale. If a
#' numeric vector of length 2, will be used as lower and upper quantiles 
#' beyond which to trim.
#' @param colors The heatmap colors to use.
#' @param scale_title The title of the scale.
#' @param title_size The size of the heatmap titles.
#' @param use_raster Passed to \code{\link[EnrichedHeatmap]{EnrichedHeatmap}}.
#' @param row_order Optional order of the rows.
#' @param cluster_rows Whether to cluster rows.
#' @param ... Passed to \code{\link[EnrichedHeatmap]{EnrichedHeatmap}}
#' @param top_annotation Either a logical indicating whether or not to plot the
#' summary profile at the top of each heatmap, or any other value passed to
#' \code{\link[EnrichedHeatmap]{EnrichedHeatmap}}.
#' @param axis_name A vector of length 3 giving the labels to put respectively
#' on the left, center and right of the x axis of each heatmap.
#'
#' @return A HeatmapList object
#' @import EnrichedHeatmap
#' @importFrom viridisLite inferno
#' @importFrom circlize colorRamp2
#' @importFrom matrixStats colSds
#' @importFrom grid gpar
#' @export
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
plotEnrichedHeatmaps <- function(ml, trim=c(0.01,0.99), colors=inferno(100),
                         scale_title="density", title_size=11, use_raster=NULL, 
                         row_order=NULL, cluster_rows=FALSE, axis_name=NULL, 
                         scale_rows=FALSE, top_annotation=TRUE, minRowVal=0, ...){
  ml <- .comparableMatrices(ml)
  stopifnot(length(trim) %in% 1:2 && all(trim>=0 & trim <=1))
  stopifnot(length(scale_rows)==1)
  if(isTRUE(minRowVal>0)){
    rMax <- do.call(cbind, lapply(ml, FUN=matrixStats::rowMaxs))
    w <- which(rMax>minRowVal)
    ml <- lapply(ml, FUN=function(x) x[w,seq_len(ncol(x))])
  }
  if(!isFALSE(scale_rows)){
    if(isTRUE(scale_rows) || scale_rows=="global"){
      m <- do.call(cbind, ml)
      rM <- rowMeans(m, na.rm=TRUE)
      rSD <- matrixStats::rowSds(m, na.rm=TRUE)
      ml <- lapply(ml, FUN=function(x) (x-rM)/rSD)
      ################
    }else{
      ml <- lapply(ml, FUN=function(x) t(scale(t(x))))
    }
  }
  if(is.null(use_raster)) use_raster <- nrow(ml[[1]])>1000
  hl <- NULL
  ylim <- c(0,max(unlist(lapply(ml,FUN=function(x){
    max(colMeans(x)+matrixStats::colSds(x)/sqrt(nrow(x)))
  }))))
  if(length(trim)==1) trim <- c(0,trim)
  trim <- sort(trim)
  common_min <- min(unlist(lapply(ml,prob=trim[1],FUN=quantile)))
  common_max <- max(unlist(lapply(ml,prob=trim[2],FUN=quantile)))
  breaks <- seq(from=common_min, to=common_max, length.out=length(colors))
  col_fun <- circlize::colorRamp2(breaks, colors)
  if(isTRUE(cluster_rows)){
    cluster_rows <- hclust(dist(do.call(cbind, ml)))
  }else if(is.null(row_order)){
    row_order <- order(-rowMeans(do.call(cbind, lapply(ml, enriched_score))))
  }
  if(is.null(axis_name) || length(axis_name)==1){
    a <- attributes(ml[[1]])$extend
    if(all((a %% 1000)==0)){
      a <- paste0(c("-","+"),a/1000,"kb")
    }else{
      a <- paste0(c("-","+"),a,"bp")
    }
    axis_name <- c(a[1], ifelse(is.null(axis_name),"center",axis_name), a[2])
  }
  for(m in names(ml)){
    isLast <- m==rev(names(ml))[1]
    if(isFALSE(top_annotation)){
      top_annotation <- NULL
    }else if(isTRUE(top_annotation)){
      top_annotation <- HeatmapAnnotation(
        enriched=anno_enriched(ylim=ylim, show_error=TRUE, axis=isLast))
    }
    hl <- hl + EnrichedHeatmap(ml[[m]], column_title=m, col=col_fun, ...,
                column_title_gp=gpar(fontsize=title_size), axis_name=axis_name,
                cluster_rows = cluster_rows, row_order=row_order,
                top_annotation=top_annotation, show_heatmap_legend=isLast,
                use_raster=use_raster, name=ifelse(isLast,scale_title,m) )
  }
  hl
}


#' meltSignals
#'
#' Aggregates and melts a list of signal matrices, for plotting (with ggplot).
#'
#' @param ml A named list of matrices as produced by \code{\link{signal2Matrix}}
#' @param fun An optional custom aggregation function (or named list thereof).
#'
#' @return A data.frame.
#' @export
#' @importFrom matrixStats colMedians
meltSignals <- function(ml, fun=NULL){
  stopifnot(is.list(ml))
  ml <- .comparableMatrices(ml, checkAttributes=TRUE)
  a <- attributes(ml[[1]])
  x <- c( -1*(a$extend[1]/length(a$upstream_index))*
            rev(seq_along(a$upstream_index)),
          (a$extend[2]/length(a$downstream_index))*
            seq_along(a$downstream_index) )
  d <- data.frame( row.names=NULL, position=rep(x, length(ml)),
                   sample=rep(factor(names(ml), names(ml)), each=length(x)) )
  d$mean <- unlist(lapply(ml, colMeans))
  d$SD <- unlist(lapply(ml, matrixStats::colSds))
  d$SE <- d$SD/sqrt(nrow(ml[[1]]))
  d$median <- unlist(lapply(ml, matrixStats::colMedians))
  if(!is.null(fun)){
    if(!is.list(fun)) fun <- list(fn.value=fun)
    stopifnot(all(unlist(lapply(fun,is.function))))
    if(is.null(names(fun)))
      names(fun) <- make.unique(rep("fn.value",length(fun)))
    for(fn in names(fun))
      d[[fn]] <- unlist(lapply(ml, FUN=function(x) apply(x,2,fun[[fn]])))
  }
  d
}

#' mergeSignalMatrices
#'
#' @param ml A named list of matrices as produced by \code{\link{signal2Matrix}}
#' @param aggregation The method to aggregate matrices
#'
#' @return A `normalizedMatrix` object
#' @export
mergeSignalMatrices <- function(ml, aggregation=c("mean","sum","median")){
  aggregation <- match.arg(aggregation)
  ml <- .comparableMatrices(ml, checkAttributes=TRUE)
  switch(aggregation,
    mean=Reduce("+", ml)/length(ml),
    sum=Reduce("+", ml),
    median=.medianSignalMat(ml)
  )
}

#' @importFrom matrixStats rowMedians
.medianSignalMat <- function(ml){
  x <- matrixStats::rowMedians(do.call(cbind,lapply(ml, as.numeric)))
  y <- ml[[1]]
  y[seq_len(nrow(y)),seq_len(ncol(y))] <- x
  y
}

#' renormalizeBorders
#'
#' This function renormalizes a list of signal matrices on the assumption that
#' the left/right borders of the matrices represent background signal which 
#' should be equal across samples.
#' \strong{This is not a safe normalization procedure}: it will work only if 
#' 1) the left/right borders of the matrices are sufficiently far from the 
#' signal (e.g. peaks), and 2) the signal-to-noise ratio is comparable across 
#' samples.
#'
#' @param ml A named matrix list as produced by \code{\link{signal2Matrix}}.
#' @param method Normalization method, passed to 
#' \code{\link[edgeR]{calcNormFactors}}.
#'
#' @return A renormalized list of signal matrices.
#' @export
#' @importFrom edgeR calcNormFactors
renormalizeBorders <- function(ml, method="TMM"){
  ml <- .comparableMatrices(ml, checkAttributes=TRUE)
  b <- do.call(cbind, lapply(ml, FUN=function(x) c(x[,1],x[,ncol(x)])))
  nf <- calcNormFactors(b, method=method, lib.size=rep(1,ncol(b)))
  ml <- rescaleSignalMatrices(ml, 1/nf)
  ml
}

#' rescaleSignalMatrices
#'
#' @param ml A named matrix list as produced by \code{\link{signal2Matrix}}.
#' @param scaleFactors A numeric vector of same length as `ml`, 
#' indicating the scaling factors by which to multiply each matrix.
#'
#' @return A renormalized list of signal matrices.
#' @export
rescaleSignalMatrices <- function(ml, scaleFactors){
  stopifnot(is.numeric(scaleFactors) && length(scaleFactors)==length(ml))
  ml <- .comparableMatrices(ml, checkAttributes=TRUE)
  for(i in seq_along(scaleFactors)) ml[[i]] <- ml[[i]]*scaleFactors[i]
  ml
}