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
#' @param column_title_gp Graphic parameters of the column titles (see 
#'   \code{\link[grid]{gpar}})
#' @param row_order Optional order of the rows.
#' @param row_split Splitting of rows.
#' @param cluster_rows Whether to cluster rows.
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
#' @param use_raster Logical; whether to render the heatmap body as a raster 
#'   image. Turned on by default if any of the matrix dimensions is greater than
#'   2000.
#' 
#' @details 
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
#' @importFrom SummarizedExperiment rowRanges
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
plotEnrichedHeatmaps <- function(ml, trim=c(0.02,0.98), assay=1L, 
                                 colors=inferno(100), scale_title="density",
                                 column_title_gp=gpar(fontsize=11), 
                                 row_order=NULL, cluster_rows=FALSE, 
                                 row_split=NULL, axis_name=NULL, minRowVal=0, 
                                 scale_rows=FALSE, top_annotation=TRUE, 
                                 left_annotation=NULL, right_annotation=NULL,
                                 mean_color="red", mean_scale_side=NULL, 
                                 show_heatmap_legend=TRUE, use_raster=NULL, ...){
  hasMean <- isTRUE(top_annotation) || is(top_annotation, "list") || 
    (is.character(top_annotation) && "mean" %in% top_annotation)
  meanPars <- list()
  if(hasMean && is(top_annotation, "list")){
    meanPars <- top_annotation
    top_annotation <- NULL
  }
  if(is(ml, "SummarizedExperiment")){
    left_annotation <- .parseRowAnn(ml, left_annotation)
    right_annotation <- .parseRowAnn(ml, right_annotation)
    top_annotation <- .parseTopAnn(ml, top_annotation)
    ml <- .ese2ml(ml, assay=assay)
  }
  ml <- .comparableMatrices(ml)
  if(is.null(use_raster))
    use_raster <- any(unlist(lapply(ml, dim))>2000)
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
  ymin <- min(c(0,unlist(lapply(ml,FUN=min))))
  ymax <- max(unlist(lapply(ml, FUN=function(x){
    max(unlist(lapply(split(seq_along(rs), rs), FUN=function(i){
      x2 <- unclass(x)[i,,drop=FALSE]
      max(colMeans(x2)+matrixStats::colSds(x2)/sqrt(nrow(x2)))
    })))
  })))
  col_fun <- .getColFun(ml, trim, colors=colors)
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
      e <- paste0(c("-","+"),round(e/1000),"kb")
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
  
  hasAnno <- !(is.null(top_annotation) && is.null(left_annotation) &&
                 is.null(right_annotation) && is.null(row_split))
  if(is.null(mean_scale_side)) mean_scale_side <- ifelse(hasAnno,"left","both")
  hl <- NULL
  for(i in seq_along(ml)){
    m <- names(ml)[i]
    isFirst <- i==1
    isLast <- m==rev(names(ml))[1]
    TA <- NULL
    if(hasMean){
      TA <- .prepAnnoEnrich(meanPars, isFirst=isFirst, isLast, gp=mean_gp,
                            ylim=c(ymin, ymax), mean_scale_side=mean_scale_side)
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
    hl <- hl + EnrichedHeatmap( ml[[m]],
     col=col_fun, ..., column_title=m, column_title_gp=column_title_gp,
     cluster_rows=cluster_rows, row_order=row_order, row_split=row_split,
     show_heatmap_legend=isLast && show_heatmap_legend,
     left_annotation=la, right_annotation=ra,
     top_annotation=TA, axis_name=axis_name, use_raste=use_raster,
     name=ifelse(isLast,scale_title,paste(scale_title,m)) )
  }
  hl
}

.getColFun <- function(ml, trim, colors){
  stopifnot(is.numeric(trim) & all(trim>=0) & all(trim<=1))
  stopifnot(length(colors)>1)
  if(length(trim)==1) trim <- c(0,trim)
  trim <- sort(trim)
  r <- c( min(unlist(lapply(ml,prob=trim[1],na.rm=TRUE,FUN=quantile))),
          max(unlist(lapply(ml,prob=trim[2],na.rm=TRUE,FUN=quantile))) )
  if(r[[1]]==r[[2]]) r <- range(unlist(lapply(ml,range)))
  if(r[[1]]==r[[2]]){
    warning("There appears to be no signal in the data!")
    r[[2]] <- r[[1]]+1
  }
  breaks <- seq(from=r[[1]], to=r[[2]], length.out=length(colors))
  circlize::colorRamp2(breaks, colors)
}

.prepAnnoEnrich <- function(par, isFirst, isLast, ...,
                            mean_scale_side=c("both","left","right","none")){
  if(is.null(par) || isFALSE(par)) return(NULL)
  if(is(par, "HeatmapAnnotation")) return(par)
  if(isTRUE(par)) par <- list()
  stopifnot(is.list(par))
  mean_scale_side <- match.arg(mean_scale_side)
  side <- "right"
  if( (!isFirst && !isLast) || mean_scale_side=="none"){
    axis <- FALSE
  }else{
    axis <- (mean_scale_side!="left" && isLast) ||
              (mean_scale_side!="right" && isFirst)
    if(isFirst) side <- "left"
  }
  defPars <- list(show_error=TRUE, axis=axis, axis_param=list(side=side), ...)
  c(defPars[setdiff(names(defPars), names(par))], par)
}

.parseRowAnn <- function(RD, fields){
  if(is(RD, "RangedSummarizedExperiment")){
    RD <- as.data.frame(rowRanges(RD))
  }else if(is(RD, "SummarizedExperiment")){
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
  if(isFALSE(fields) || isTRUE(fields)) return(NULL)
  if(is.list(CD)) return(NULL)
  if(is(CD, "SummarizedExperiment")) CD <- colData(CD)
  fields <- intersect(fields, colnames(CD))
  if(length(fields)==0) return(NULL)
  as.data.frame(CD[,fields,drop=FALSE])
}

.prepTopAnn <- function(df, i, ml){
  print(df)
  if(is.null(df)) return(NULL)
  df <- df[rep(i,ncol(ml[[i]])),,drop=FALSE]
  row.names(df) <- NULL
  df
}