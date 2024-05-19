#' plotEnrichedHeatmaps: Plots heatmaps of signals around a set of regions
#' 
#' Plots enrichment heatmaps from the output of \code{\link{signal2Matrix}} 
#' (i.e. an EnrichmentSE object or a list of signal matrices). This is a 
#' convenience wrapper around \code{\link[EnrichedHeatmap]{EnrichedHeatmap}}.
#'
#' @param ml A named matrix list as produced by \code{\link{signal2Matrix}}, or
#'  an `EnrichmentSE` object.
#' @param trim The quantile above which to trim values for the colorscale. If a
#' numeric vector of length 2, will be used as lower and upper quantiles 
#' beyond which to trim.
#' @param colors The heatmap colors to use, a vector of at least two colors 
#'   between which to interpolate. Can also be a list of such color scales, with
#'   as many slots as there are tracks in `ml`. If a list of single colors, a
#'   color scale from white to that color will be used for each track. Defaults
#'   to the 'inferno' viridis palette.
#' @param multiScale Logical; whether to use a different scale for each track.
#'   Defaults to TRUE is `colors` is a list, otherwise FALSE.
#' @param scale_title The title of the scale. Ignored if `multiScale=TRUE`.
#' @param column_title_gp Graphic parameters of the column titles (see 
#'   \code{\link[grid]{gpar}})
#' @param row_order Optional order of the rows.
#' @param row_split Variable according to which the rows should be split. This 
#'   should either be the name of a column of `rowData(ml)`, or a factor/
#'   character vector of length equal to the number of regions in `ml`.
#' @param cluster_rows Logical; whether to cluster rows. Alternatively, 
#'   `cluster_rows="sort"` will sort rows using the angle on an MDS based on the
#'   \code{\link[EnrichedHeatmap]{enriched_score}} of the signals (can be very
#'   long to compute on large matrices). `cluster_rows=FALSE` (default) results 
#'   in the traditional sorting by decreasing `enriched_score`.
#' @param ... Passed to \code{\link[EnrichedHeatmap]{EnrichedHeatmap}}
#' @param top_annotation Either a logical indicating whether or not to plot the
#' summary profile at the top of each heatmap, a named list of parameters to be 
#'   passed to `anno_enrich`, or a 
#'   \code{\link[ComplexHeatmap]{HeatmapAnnotation-class}} object that will be 
#'   passed to \code{\link[EnrichedHeatmap]{EnrichedHeatmap}}. Additionally, if
#'   `ml` is a `ESE` object, `top_annotation` can be a vector of colData column
#'   names.
#' @param axis_name A vector of length 3 giving the labels to put respectively
#' on the left, center and right of the x axis of each heatmap.
#' @param assay Assay to use (ignored unless `ml` is an ESE object)
#' @param minRowVal Minimum value a row should have to be included
#' @param scale_rows Whether to scale rows, either FALSE (default), 'local' 
#'   (scales each matrix separately) or 'global'.
#' @param left_annotation Passed to \code{\link[EnrichedHeatmap]{EnrichedHeatmap}}
#' @param right_annotation Passed to \code{\link[EnrichedHeatmap]{EnrichedHeatmap}}
#' @param show_heatmap_legend Logical, whether to show the heatmap legend
#' @param mean_color Color of the mean signal line in the top annotation. If
#'   `row_split` is used, `mean_color` can be a named vector indicating the 
#'   colors for each cluster. Can also be a `gpar` object.
#' @param mean_scale_side The side on which to show the y-axis scale of the mean
#'   plots. Either "both" (default), "left", "right", or "none".
#' @param mean_trim Logical; whether to apply the trimming also to the mean plot.
#' @param use_raster Logical; whether to render the heatmap body as a raster 
#'   image. Turned on by default if any of the matrix dimensions is greater than
#'   1500.
#' 
#' @details 
#' When plotting large matrices, the heatmap body will be rasterized to keep its
#' memory footprint decent. For this to work well, make sure the `magick` 
#' package is installed. Depending on your settings, if the heatmap is very
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
#' @importFrom SummarizedExperiment rowRanges
#' @export
#' @examples 
#' # we first load an example EnrichmentSE, as produced by signal2Matrix:
#' data(exampleESE)
#' plotEnrichedHeatmaps(exampleESE)
#' # we could also just plot one with:
#' # plotEnrichedHeatmaps(exampleESEm[,1])
#' # or change the aesthetics, e.g.:
#' plotEnrichedHeatmaps(exampleESE, trim=0.98, scale_title="RPKM", 
#'                      colors=c("white","darkred"), row_title="My regions")
#' # any argument accepted by `EnrichedHeatmap` (and hence by 
#' # `ComplexHeatmap::Heatmap`) can be used.
plotEnrichedHeatmaps <- function(ml, trim=c(0.02,0.98), assay=1L, 
                                 colors=NULL, scale_title="density",
                                 column_title=NULL, multiScale=NULL,
                                 column_title_gp=gpar(fontsize=11), 
                                 row_order=NULL, cluster_rows=FALSE, 
                                 row_split=NULL, axis_name=NULL, minRowVal=0, 
                                 scale_rows=FALSE, top_annotation=TRUE, 
                                 left_annotation=NULL, right_annotation=NULL,
                                 mean_color="red", mean_scale_side=NULL, 
                                 mean_trim=TRUE, show_heatmap_legend=TRUE,
                                 use_raster=NULL, ...){
  hasMean <- isTRUE(top_annotation) || is(top_annotation, "list") || 
    (is.character(top_annotation) && "mean" %in% top_annotation)
  # avoid inconsistencies when resizing the matrix
  if(hasMean & isFALSE(mean_trim) && 
     !is.null(rrm <- list(...)[["raster_resize_mat"]]) && rrm){
    warning("Enforcing mean_trim=TRUE due to the usage of 'raster_resize_mat'")
    mean_trim <- TRUE
  }
  meanPars <- list()
  if(is(top_annotation, "list")){
    meanPars <- top_annotation
    top_annotation <- NULL
  }else if(is.logical(top_annotation)){
    top_annotation <- NULL
  }
  if(is(ml, "EnrichmentSE")){
    if(is.null(colors) && !isFALSE(multiScale) && !is.null(ml$hmcolors) &&
       is.list((ml$hmcolors))) colors <- ml$hmcolors
    left_annotation <- .parseRowAnn(ml, left_annotation)
    right_annotation <- .parseRowAnn(ml, right_annotation)
    top_annotation <- .parseTopAnn(ml, top_annotation)
    if(!is.null(row_split) && is.character(row_split) && length(row_split)==1){
      if(!(row_split %in% colnames(rowData(ml))))
        stop("`row_split` value not found in the rowData columns.")
      row_split <- rowData(ml)[,row_split]
    }
    ml <- .ese2ml(ml, assay=assay)
  }
  if(is.null(colors)) colors <- inferno(100)
  ml <- .comparableMatrices(ml)
  if(is.null(use_raster))
    use_raster <- any(unlist(lapply(ml, dim))>1500)
  stopifnot(length(trim) %in% 1:2 && all(trim>=0 & trim <=1))
  stopifnot(length(scale_rows)==1)
  if(!is.null(row_split)) stopifnot(length(row_split)==nrow(ml[[1]]))
  if(isTRUE(minRowVal>0)){
    rMax <- matrixStats::rowMaxs(do.call(cbind,
                                         lapply(ml, FUN=matrixStats::rowMaxs)))
    w <- which(rMax>minRowVal)
    ml <- lapply(ml, i=w, FUN=.resizeNmatrix)
    if(!is.null(row_split)) row_split <- row_split[w]
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
  mean_gp <- gpar(col="red")
  if(is.null(row_split)){
    rs <- rep(1L, nrow(ml[[1]]))
  }else{
    rs <- row_split <- as.factor(row_split)
    if(is(mean_color, "gpar")){
      mean_gp <- mean_color
    }else if(length(mean_color)>1){
      if(is.null(names(mean_color)))
        stop("When clustering rows, provide the names of the clusters in `mean_color`.")
      stopifnot(all(levels(row_split) %in% names(mean_color)))
      mean_gp <- gpar(col=mean_color[levels(row_split)])
    }
  }
  
  if(is.null(multiScale)) multiScale <- is.list(colors)
  if(is.list(colors) && length(colors)!=length(ml))
      stop("The number of input color scales should either be one or the ",
           "number of tracks.")
  if(multiScale){
    if(!is.list(colors)){
      colors <- lapply(seq_along(ml), FUN=function(i) colors)
    }else{
      colors <- lapply(colors, FUN=function(x){
        if(length(x)==1) x <- c("white",x)
        x
      })
    }
    ml_trimmed <- lapply(ml, trim=trim, FUN=.applyTrimming)
    col_fun <- lapply(seq_along(ml), FUN=function(i){
      .getColFun(ml[[i]], trim, colors=colors[[i]])
    })
  }else{
    ml_trimmed <- .applyTrimming(ml, trim)
    col_fun <- .getColFun(ml, trim, colors=colors)
    col_fun <- lapply(seq_along(ml), FUN=function(i) col_fun)
  }
  if(isTRUE(mean_trim)){
    ml <- ml_trimmed # don't just trim the scale
    trim <- 1
  }

  ymin <- min(c(0,unlist(lapply(ml, na.rm=TRUE, FUN=min))))
  ymax <- max(unlist(lapply(ml, FUN=function(x){
    max(unlist(lapply(split(seq_along(rs), rs), FUN=function(i){
      x2 <- unclass(x)[i,,drop=FALSE]
      x2[is.na(x2)] <- 0
      max(colMeans(x2)+matrixStats::colSds(x2)/sqrt(nrow(x2)))
    })))
  })))
  if(isTRUE(cluster_rows)){
    cluster_rows <- hclust(dist(do.call(cbind, ml_trimmed)))
  }else if(is.null(row_order)){
    es <- do.call(cbind, lapply(ml_trimmed, enriched_score))
    if(cluster_rows=="sort"){
      row_order <- .mdsSortRows(es)
    }else{
      row_order <- order(-rowMeans(es))
    }
  }
  if(!is.null(row_order)) cluster_rows <- FALSE
  
  hasAnno <- !(is.null(top_annotation) && is.null(left_annotation) &&
                 is.null(right_annotation) && is.null(row_split))
  if(is.null(mean_scale_side)) mean_scale_side <- ifelse(hasAnno,"left","both")
  hl <- NULL
  for(i in seq_along(ml)){
    m <- names(ml)[i]
    axisLabels <- .getAxisLabels(ml[[i]], axis_name=axis_name)
    isFirst <- i==1
    isLast <- m==rev(names(ml))[1]
    TA <- NULL
    if(hasMean){
      if(multiScale){
        TA <- .prepAnnoEnrich(meanPars, isFirst=TRUE, FALSE, gp=mean_gp,
                              axis=FALSE)
      }else{
        TA <- .prepAnnoEnrich(meanPars, isFirst=isFirst, isLast, gp=mean_gp,
                              ylim=c(ymin, ymax), mean_scale_side=mean_scale_side)
      }
      if(is.null(top_annotation)){
        TA <- HeatmapAnnotation(enriched=do.call(anno_enriched, TA))
      }else{
        pars <- c( as.list(.prepTopAnn(top_annotation,i,ml)),
                   list(
                     enriched=do.call(anno_enriched, TA),
                     show_annotation_name=(isFirst || isLast),
                     annotation_name_side=ifelse(isLast,"right","left")
                   ))
        TA <- do.call(HeatmapAnnotation, pars)
      }
    }else if(!is.null(top_annotation)){
      TA <- HeatmapAnnotation(df=.prepTopAnn(top_annotation,i,ml),
                              show_annotation_name=(isFirst || isLast),
                              annotation_name_side=ifelse(isLast,"right","left"))
    }
    la <- NULL
    if(isFirst){
      if(is.null(la <- left_annotation) &&
       !is.null(row_split) && !is(mean_color, "gpar") && length(mean_color)>1){
        la <- rowAnnotation(cluster=row_split, annotation_name_side="bottom",
                            col=list(cluster=mean_color), show_legend=FALSE)
      }else if(is.data.frame(la)){
        la <- rowAnnotation(df=la, annotation_name_side="bottom")
      }
    }
    ra <- NULL
    if(isLast && !is.null(ra <- right_annotation)){
      if(is.data.frame(ra))
        ra <- rowAnnotation(df=ra, annotation_name_side="bottom")
    }
    if(multiScale){
      hmname <- m
    }else{
      hmname <- ifelse(isLast,scale_title,paste(scale_title,m,sep="\n"))
    }
    clt <- ifelse(is.null(column_title),m,column_title)
    hl <- hl + EnrichedHeatmap( ml[[m]], ..., col=col_fun[[i]], 
       column_title=clt, column_title_gp=column_title_gp,
       cluster_rows=cluster_rows, row_order=row_order, row_split=row_split,
       show_heatmap_legend=(multiScale || isLast) && show_heatmap_legend,
       left_annotation=la, right_annotation=ra,
       top_annotation=TA, axis_name=axisLabels, use_raste=use_raster,
       name=hmname )
  }
  hl
}

.getAxisLabels <- function(m, axis_name=NULL){
  a <- attributes(m)
  neededAxisLabs <- ifelse(length(a$target_index)==0,3,4)
  if(is.null(axis_name) || length(axis_name)<neededAxisLabs){
    e <- paste0(c("-","+"), formatGenomicDist(a$extend))
    if(length(a$target_index)>0){
      if(is.null(axis_name)) axis_name <- c("start","end")
      axis_name <- head(axis_name,2)
      axis_name <- c(axis_name,c("start","end")[-seq_along(axis_name)])
    }else{
      axis_name <- ifelse(is.null(axis_name),"center",axis_name)
    }
    axis_name <- c(e[1], axis_name, e[2])
  }
  axis_name
}



.getColFun <- function(ml, trim, colors){
  stopifnot(length(colors)>1)
  r <- .getTrimPoints(ml, trim)
  breaks <- seq(from=r[[1]], to=r[[2]], length.out=length(colors))
  circlize::colorRamp2(breaks, colors)
}

# get range beyond which to trim
.getTrimPoints <- function(ml, trimQ){
  stopifnot(is.numeric(trimQ) & all(trimQ>=0) & all(trimQ<=1))
  if(length(trimQ)==1) trimQ <- c(0,trimQ)
  trimQ <- sort(trimQ)
  if(!is.list(ml)) ml <- list(ml)
  tryCatch({
    r <- c( min(unlist(lapply(ml,prob=trimQ[1],na.rm=TRUE,FUN=quantile))),
            max(unlist(lapply(ml,prob=trimQ[2],na.rm=TRUE,FUN=quantile))) )
    if(r[[1]]==r[[2]]) r <- range(unlist(lapply(ml,range)))
    if(r[[1]]==r[[2]]){
      warning("There appears to be no signal in the data!")
      r[[2]] <- r[[1]]+1
    }
    r
  }, error=function(e){
    stop("Unable to apply trimming - this typically happens when the input ",
         "contains less non-zero values than the trimmed range.\nDouble-check ",
         "your object or try to reduce the trimming.")
  })
  r
}

# trim values beyond trim range in signal matrices
.applyTrimming <- function(ml, trim){
  if(singleM <- !is.list(ml)){
    ml <- list(ml)
  }
  br <- .getTrimPoints(ml, trim)
  ml <- lapply(ml, FUN=function(x){
    x[which(x>br[2])] <- br[2]
    x[which(x<br[1])] <- br[1]
    x
  })
  if(singleM) ml <- ml[[1]]
  ml
}

.prepAnnoEnrich <- function(par, isFirst, isLast, axis=NULL, ...,
                            mean_scale_side=c("both","left","right","none")){
  if(is.null(par) || isFALSE(par)) return(NULL)
  if(is(par, "HeatmapAnnotation")) return(par)
  if(isTRUE(par)) par <- list()
  stopifnot(is.list(par))
  mean_scale_side <- match.arg(mean_scale_side)
  side <- "right"
  if(is.null(axis)){
    if( (!isFirst && !isLast) || mean_scale_side=="none"){
      axis <- FALSE
    }else{
      axis <- (mean_scale_side!="left" && isLast) ||
                (mean_scale_side!="right" && isFirst)
      if(isFirst) side <- "left"
    }
  }
  defPars <- list(show_error=TRUE, axis=axis, axis_param=list(side=side), ...)
  c(defPars[setdiff(names(defPars), names(par))], par)
}

.parseRowAnn <- function(RD, fields){
  if(inherits(RD, "RangedSummarizedExperiment")){
    RD <- as.data.frame(rowRanges(RD))
  }else if(inherits(RD, "SummarizedExperiment")){
    RD <- rowData(RD)
  }
  if(!is.character(fields) || length(fields)==nrow(RD)) return(fields)
  if(length(fields)==0) return(NULL)
  fields <- intersect(fields, colnames(RD))
  if(length(fields)==0){
    warning("Unknown fields in row annotation!")
    return(NULL)
  }
  as.data.frame(RD[,fields,drop=FALSE])
}

.parseTopAnn <- function(CD, fields){
  if(is.null(fields) || isFALSE(fields) || isTRUE(fields)) return(NULL)
  if(is.list(CD)) return(NULL)
  if(inherits(CD, "SummarizedExperiment")) CD <- colData(CD)
  fields <- intersect(fields, colnames(CD))
  if(length(fields)==0) return(NULL)
  as.data.frame(CD[,fields,drop=FALSE])
}

.prepTopAnn <- function(df, i, ml){
  if(is.data.frame(df) || is(df,"DFrame")){
    df <- df[rep(i,ncol(ml[[i]])),,drop=FALSE]
    row.names(df) <- NULL
  }
  df
}
