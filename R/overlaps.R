#' regionsToUpset
#' 
#' Prepares sets of regions for UpSet overlap representation.
#' A wrapper around \code{\link[UpSetR]{upset}} for comparing multiple sets of
#' genomic ranges.
#'
#' @param x A named list of genomic ranges (or paths to bed files)
#' @param reference The method for creating the reference windows ('reduce' or
#'   'disjoin'). Alternatively, a `GRanges` object of reference windows.
#' @param returnList Logical; whether to return the list of regions instead of 
#'   plotting.
#' @param ... Further arguments specifying how the overlaps are done, 
#' passed to \code{\link[GenomicRanges]{findOverlaps-methods}}).
#'
#' @return A data.frame of set inclusions which can be directly input to
#'   \code{\link[ComplexHeatmap]{make_comb_mat}}, and then
#'   \code{\link[ComplexHeatmap]{UpSet}}.
#' @export
#' @importFrom GenomicRanges reduce disjoin 
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom ComplexHeatmap UpSet make_comb_mat
#'
#' @examples
#' # random list of GRanges:
#' grl <- lapply(c(A=10,B=20,C=30), FUN=function(x){
#'   GRanges("seq1", IRanges(runif(x,1,1000), width=20))
#' })
#' input_for_upset <- regionsToUpset(grl)
#' # we would then plot the data with:
#' ComplexHeatmap(UpSet(make_comb_mat(input_for_upset)))
regionsToUpset <- function(x, reference=c("reduce","disjoin"), returnList=FALSE,
                           ignore.strand=FALSE, maxgap=-1L, minoverlap=0L, ...){
  if(is.character(x)) x <- as.list(x)
  stopifnot(length(x)>1)
  if(is.list(x)){
    if(is.null(names(x)) && all(unlist(lapply(x,is.character))) && 
       all(lengths(x)==1))
      names(x) <- unlist(x)
    if(any(unlist(lapply(x,FUN=function(x){
      if(is.character(x)){
        if(length(x)!=0) return(TRUE)
        if(!file.exists(x)) return(TRUE)
      }else if(!is(x,"GRanges")){
        return(TRUE)
      }
      FALSE
    })))){
      stop("`x` should either be i) a named list of GRanges, ii) a 
            GRangesList, or iii) a named list of character vectors, each of
            length 1 and providing the path to a bed file.")
    }
    x <- lapply(x, FUN=function(x){
      if(is.character(x)){
        if(grepl("\\.rds$",x,ignore.case=TRUE)){
          y <- readRDS(x)
        }else{
          y <- importBedlike(x)
        }
        if(!is(y,"GRanges"))
          stop(y," does not appear to contain genomic ranges.")
        x <- y
      }
      x
    })
    if(is.list(x) && all(unlist(lapply(x,class2="GRanges",is))))
      stopifnot(all(unlist(lapply(x,class2="GRanges",is))))
    x <- as(x, "GRangesList")
  }
  stopifnot(is(x,"GRangesList"))
  if(!is(reference,"GRanges"))
    reference <- switch(match.arg(reference),
                        reduce=reduce(unlist(x),ignore.strand=ignore.strand),
                        disjoin=disjoin(unlist(x),ignore.strand=ignore.strand))
  x <- vapply(x, FUN.VALUE=logical(length(reference)), FUN=function(x){
    overlapsAny(reference, x, ignore.strand=ignore.strand, 
                maxgap=maxgap, minoverlap=minoverlap)
  })
  if(returnList) return(apply(x, 2, function(x) reference[which(x)]))
  as.data.frame(row.names=as.character(reference), x)
}


#' regionOverlaps
#'
#' A wrapper for visualizing pairwise-wise overlaps across multiple sets of
#' genomic ranges.
#' 
#' @param listOfRegions A named list of two or more (non-empty) `GRanges`
#' @param mode Either 'reduced' or 'pairwise'. 'reduced' first uses `reduce` to
#'   get a set of reference regions which are, based on overlap, contained or 
#'   not in the different sets. It is thus symmetrical. `pairwise` does pairwise
#'   overlap between the sets of regions; it is asymmetrical and slower to 
#'   compute.
#' @param colorBy Whether to color by 'overlapCoef' (default), or by 'jaccard'
#'   index.
#' @param ignore.strand Logical; whether to ignore strand for overlaps (default
#'   TRUE).
#' @param color Heatmap colorscale
#' @param cluster Logical; whether to cluster rows/columns
#' @param number_color Values color
#' @param ... Passed to \code{\link[ComplexHeatmap]{pheatmap}}
#'
#' @return A `Heatmap` showing the overlap coefficient as colors, and the 
#'   overlap size as values.
#' @importFrom IRanges overlapsAny
#' @importFrom ComplexHeatmap pheatmap
#' @importFrom viridisLite plasma
#' @export
#' @examples
#' # random list of GRanges:
#' grl <- lapply(c(A=10,B=20,C=30), FUN=function(x){
#'   GRanges("seq1", IRanges(runif(x,1,1000), width=20))
#' })
#' regionOverlaps(grl)
regionOverlaps <- function(listOfRegions, mode=c("reduced","pairwise"),
                           ignore.strand=TRUE, cluster=length(listOfRegions)>2,
                           colorBy=c("overlapCoef","jaccard"), ...,
                           color=viridis::plasma(100), number_color="black"){
  stopifnot(length(listOfRegions)>1 && all(lengths(listOfRegions)>0) &&
              all(sapply(listOfRegions,class2="GRanges",is)))
  mode <- match.arg(mode)
  colorBy <- match.arg(colorBy)
  if(mode=="reduced"){
    r <- reduce(GRangesList(listOfRegions), ignore.strand=ignore.strand)
    m <- sapply(listOfRegions, \(x) overlapsAny(r, m,
                                                ignore.strand=ignore.strand))
    o <- t(m) %*% m
  }else{
    o <- suppressWarnings(sapply(listOfRegions, FUN=function(x){
      sapply(listOfRegions, FUN=function(y){
        if(identical(x,y)) return(length(x))
        sum(overlapsAny(x,y,ignore.strand=ignore.strand))
      })
    }))
  }
  sizes <- lengths(listOfRegions)
  if(colorBy=="overlapCoef"){
    co <- o/outer(sizes, sizes, "min")
  }else{
    co <- o/(outer(sizes, sizes, "+")-o)
  }
  diag(o) <- NA_real_
  h <- NULL
  if(isTRUE(cluster)) h <- hclust(dist(co))
  if(is(cluster,"hclust") || is(cluster,"dendrogram")) h <- cluster
  if(length(unique(as.numeric(co)))<3) co[is.na(co)] <- 1
  dimnames(co) <- dimnames(o)
  cn <- ifelse(colorBy=="jaccard", "Jaccard", "overlap\ncoefficient")
  ComplexHeatmap::pheatmap(co, display_numbers=o, number_color=number_color,
                           cluster_rows=h, cluster_cols=h,
                           name=cn, color=color, ...)
}



#' regionCAT
#'
#' Computes/plots the 'concordance at the top' (CAT) of two lists of genomic
#'   regions.
#'
#' @param regions1,regions2 A GRanges object with a `score` metadata column according to 
#'   which the regions will be ranked (descending).
#' @param start The rank at which to start plotting (removes large variations
#'   at the beginning when very few regions are considered)
#' @param concord.type Concordance type to plot, either 'inTop', 'inAll', or 
#'   'both' (see details). Ignored if `returnData=TRUE`.
#' @param returnData Logical; whether to return the data instead of plotting.
#' @param ignore.strand Logical; whether to ignore the strand for computing
#'   overlap (default TRUE)
#'   
#' @details 
#' The two concordance types are as follows:
#' * 'inTop' indicates the proportion of the top X regions that are in the top 
#'   X in both lists.
#' * 'all' indicates the proportion of the top X regions that are anywhere in 
#'   the other list (since this relationship is asymmetrical, the mean of both
#'   two directions is used).
#'
#' @return A ggplot object, or a data.frame if `returnData=TRUE`.
#' @export
#' @importFrom GenomicRanges reduce
#' @examples
#' # we create two GRanges with scores, which have similar high-score peaks but
#' # the rest random:
#' gr1 <- GRanges("seq1", IRanges(runif(20,1,2000), width=20),
#'                score=20:1)
#' gr2 <- GRanges("seq1", c(head(ranges(gr1),5),
#'                          IRanges(runif(15,1,2000), width=20)),
#'                score=c(20:16, sample.int(15)))
#' regionCAT(gr1,gr2)
regionCAT <- function(regions1, regions2, start=5L,
                      concord.type=c("both","inTop","inAll"),
                      returnData=FALSE, ignore.strand=TRUE){
  stopifnot(is(regions1,"GRanges") && is(regions1,"GRanges"))
  stopifnot(length(regions1)>start && length(regions1)>start)
  stopifnot(!is.null(regions1$score) && !is.null(regions2$score))
  concord.type <- match.arg(concord.type)
  o <- reduce(c(regions1,regions2), with.revmap=TRUE,
              ignore.strand=ignore.strand)$revmap
  o <- unlist(o[lengths(o)>1])
  o1 <- o[o<=length(regions1)]
  o2 <- o[o>length(regions1)]-length(regions1)
  regions1$overlaps <- seq_along(regions1) %in% o1
  regions2$overlaps <- seq_along(regions2) %in% o2
  regions1 <- head(regions1$overlaps[order(-regions1$score)], 
                   min(length(regions1),length(regions2)))
  regions2 <- head(regions2$overlaps[order(-regions2$score)], 
                   min(length(regions1),length(regions2)))
  d <- data.frame(rank=seq_along(regions1), 
                  p.all=(cumsum(regions1)/seq_along(regions1) + 
                           cumsum(regions2)/seq_along(regions2))/2,
                  p.top=cumsum(regions1 & regions2)/seq_along(regions1))
  if(returnData) return(d)
  
  d <- d[d$rank>=start,]
  if(concord.type=="both"){
    d <- rbind( cbind(type=rep("inTop",nrow(d)),
                      setNames(d[,c("rank","p.top")], c("rank","prop"))),
                cbind(type=rep("inAll",nrow(d)),
                      setNames(d[,c("rank","p.all")], c("rank","prop"))) )
    return(ggplot(d, aes(rank, prop, colour=type)) + geom_line(size=1.5) +
             labs(x="Rank", y="Proportion of overlap"))
  }else if(concord.type=="inTop"){
    d$prop <- d$p.top
  }else{
    d$prop <- d$p.all
  }
  requireNamespace("ggplot2")
  ggplot2::ggplot(d, ggplot2::aes(rank, prop)) + ggplot2::geom_line(size=1.5) +
    ggplot2::labs(x="Rank", y="Proportion of overlap")
}

#' colOverlaps
#' 
#' Computes pairwise overlap metric between columns
#'
#' @param x A logical matrix.
#' @param y An optional nother logical matrix or vector. If NULL (default), 
#'   overlaps are computed between columns of `x`.
#' @param metric Either 'overlap', 'jaccard', or 'overlapCoef'.
#'
#' @returns A matrix of pairwise metric values.
#' @export
#'
#' @examples
#' m <- matrix(sample(c(TRUE,FALSE),12,replace=TRUE), nrow=4)
#' colOverlaps(m, metric="jaccard")
colOverlaps <- function(x, y=NULL, metric=c("overlap","jaccard","overlapCoef")){
  metric <- match.arg(metric)
  x <- as(x, "lMatrix")
  ys <- xs <- colSums(x)
  if(is.null(y)){
    y <- x
  }else{
    if(is.vector(y)) y <- as.matrix(y)
    y <- as(y, "lMatrix")
    ys <- colSums(y)
  }
  intersections <- t(x) %*% as.matrix(y)
  if(metric=="overlap"){
    out <- intersections
  }else if(metric=="overlapCoef"){
    minSize <- outer(xs, ys, "min")
    out <- intersections/minSize
  }else{
    unions <- outer(xs, ys, "+") - intersections
    out <- intersections/unions
  }
  if(is(out, "dgeMatrix")) out <- as.matrix(out)
  out
}
