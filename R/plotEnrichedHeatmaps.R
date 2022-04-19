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
#' @details 
#' 
#' When plotting large matrices, the heatmap body will be rasterized to keep its
#' memory footprint decent. Depending on your settings, if the heatmap is very
#' big you might hit the preset limits of `magick` base rasterization, which
#' could result in an error such as 'Image must have at least 1 frame to write a 
#' bitmap'. In such cases, you might have to degrade to a lower-quality 
#' rasterization using the additional arguments 
#' `raster_by_magick=FALSE, raster_device="CairoJPEG"` .
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
plotEnrichedHeatmaps <- function(ml, trim=c(0.02,0.98), assay=1L, colors=inferno(100),
                                 scale_title="density", title_size=11, use_raster=NULL, 
                                 row_order=NULL, cluster_rows=FALSE, axis_name=NULL, 
                                 scale_rows=FALSE, top_annotation=TRUE, minRowVal=0, 
                                 ...){
  if(is(ml, "SummarizedExperiment")) ml <- .ese2ml(ml, assay=assay)
  ml <- .comparableMatrices(ml)
  stopifnot(length(trim) %in% 1:2 && all(trim>=0 & trim <=1))
  stopifnot(length(scale_rows)==1)
  if(isTRUE(minRowVal>0)){
    rMax <- matrixStats::rowMaxs(do.call(cbind,
                                         lapply(ml, FUN=matrixStats::rowMaxs)))
    w <- which(rMax>minRowVal)
    ml <- lapply(ml, i=w, FUN=.resizeNmatrix)
  }
  if(!isFALSE(scale_rows)){
    if(isTRUE(scale_rows) || scale_rows=="global"){
      m <- do.call(cbind, ml)
      rM <- rowMeans(m, na.rm=TRUE)
      rSD <- matrixStats::rowSds(m, na.rm=TRUE)
      rSD[rSD==0] <- max(rM)
      ml <- lapply(ml, FUN=function(x) (x-rM)/rSD)
    }else{
      ml <- lapply(ml, FUN=function(x) t(scale(t(x))))
    }
  }
  if(is.null(use_raster)) use_raster <- nrow(ml[[1]])>1000
  hl <- NULL
  ymin <- min(c(0,unlist(lapply(ml,FUN=min))))
  ymax <- max(unlist(lapply(ml,FUN=function(x){
    max(colMeans(x)+matrixStats::colSds(x)/sqrt(nrow(x)))
  })))
  if(length(trim)==1) trim <- c(0,trim)
  trim <- sort(trim)
  common_min <- min(unlist(lapply(ml,prob=trim[1],na.rm=TRUE,FUN=quantile)))
  common_max <- max(unlist(lapply(ml,prob=trim[2],na.rm=TRUE,FUN=quantile)))
  breaks <- seq(from=common_min, to=common_max, length.out=length(colors))
  col_fun <- circlize::colorRamp2(breaks, colors)
  if(isTRUE(cluster_rows)){
    cluster_rows <- hclust(dist(do.call(cbind, ml)))
  }else if(is.null(row_order)){
    row_order <- order(-rowMeans(do.call(cbind, lapply(ml, enriched_score))))
  }
  a <- attributes(ml[[1]])
  neededAxisLabs <- ifelse(length(a$target_index)==0,3,4)
  if(is.null(axis_name) || length(axis_name)<neededAxisLabs){
    e <- a$extend
    if(all((e %% 1000)==0)){
      e <- paste0(c("-","+"),e/1000,"kb")
    }else{
      e <- paste0(c("-","+"),e,"bp")
    }
    if(length(a$target_index)>0){
      if(is.null(axis_name)) axis_name <- c("start","end")
      axis_name <- head(axis_name,2)
      axis_name <- c(axis_name,c("start","end")[-seq_along(axis_name)])
    }else{
      axis_name <- ifelse(is.null(axis_name),"center",axis_name)
    }
    axis_name <- c(e[1], axis_name, e[2])
  }
  for(m in names(ml)){
    isLast <- m==rev(names(ml))[1]
    if(isFALSE(top_annotation)){
      top_annotation <- NULL
    }else if(isTRUE(top_annotation)){
      top_annotation <- HeatmapAnnotation(
        enriched=anno_enriched(ylim=c(ymin,ymax), show_error=TRUE, axis=isLast))
    }
    hl <- hl + EnrichedHeatmap(ml[[m]], column_title=m, col=col_fun, ...,
                               column_title_gp=gpar(fontsize=title_size), axis_name=axis_name,
                               cluster_rows = cluster_rows, row_order=row_order,
                               top_annotation=top_annotation, show_heatmap_legend=isLast,
                               use_raster=use_raster, name=ifelse(isLast,scale_title,m) )
  }
  hl
}