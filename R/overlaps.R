#' regionUpset
#' 
#' A wrapper around \code{\link[UpSetR]{upset}} for GRanges.
#'
#' @param x A named list of genomic ranges (or paths to bed files)
#' @param reference The method for creating the reference windows ('reduce' or
#'   'disjoin'). Alternatively, a `GRanges` object of reference windows.
#' @param returnList Logical; whether to return the list instead of plotting.
#' @param ignore.strand Logical; whether to ignore strands when computing 
#' overlaps (default FALSE). Strand information is ignored if either of the 
#' compared sets of regions is unstranded.
#' @param maxgap Maximum gap between regions to count as an overlap (see 
#'  \code{\link[GenomicRanges]{findOverlaps-methods}}).
#' @param minoverlap Minimum overlap to count as a match (see 
#'  \code{\link[GenomicRanges]{findOverlaps-methods}}).
#' @param ... Further plotting arguments passed to \code{\link[UpSetR]{upset}}.
#'
#' @return A plot
#' @export
#' @importFrom UpSetR upset fromList
#' @importFrom GenomicRanges reduce disjoin 
#' @importFrom IRanges IRanges overlapsAny
#'
#' @examples
#' # random list of GRanges:
#' grl <- lapply(c(A=10,B=20,C=30), FUN=function(x){
#'   GRanges("seq1", IRanges(runif(x,1,1000), width=20))
#' })
#' regionUpset(grl)
regionUpset <- function(x, reference=c("reduce","disjoin"), returnList=FALSE,
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
          stop(paste0(y," does not appear to contain genomic ranges."))
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
  x <- lapply(x, FUN=function(x)
    which(overlapsAny(reference, x, ignore.strand=ignore.strand, 
                      maxgap=maxgap, minoverlap=minoverlap)))
  if(returnList) return(x)
  UpSetR::upset(UpSetR::fromList(x), ...)
}


#' regionOverlaps
#'
#' @param listOfRegions A named list of two or more (non-empty) `GRanges`
#' @param ignore.strand Logical; whether to ignore strand for overlaps
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
regionOverlaps <- function(listOfRegions, ignore.strand=TRUE, 
                           cluster=length(listOfRegions)>2,
                           color=viridis::plasma(100),
                           number_color="black", ...){
  stopifnot(length(listOfRegions)>1 && all(lengths(listOfRegions)>0) &&
              all(sapply(listOfRegions,class2="GRanges",is)))
  o <- suppressWarnings(sapply(listOfRegions, FUN=function(x){
    sapply(listOfRegions, FUN=function(y){
      if(identical(x,y)) return(length(x))
      sum(overlapsAny(x,y,ignore.strand=ignore.strand))
    })
  }))
  co <- sapply(seq_along(listOfRegions), FUN=function(x){
    sapply(seq_along(listOfRegions), FUN=function(y){
      if(identical(x,y)) return(NA_real_)
      round(o[x,y]/min(lengths(listOfRegions[c(x,y)])),2)
    })
  })
  h <- NULL
  if(isTRUE(cluster)) h <- hclust(dist(co))
  if(is(cluster,"hclust") || is(cluster,"dendrogram")) h <- cluster
  if(length(unique(as.numeric(co)))<3) co[is.na(co)] <- 1
  dimnames(co) <- dimnames(o)
  ComplexHeatmap::pheatmap(co, display_numbers=o, number_color=number_color,
                           cluster_rows=h, cluster_cols=h,
                           name="overlap\ncoefficient", color=color, ...)
}



#' regionCAT
#'
#' Computes/plots the 'concordance at the top' (CAT) of two lists of genomic
#'   regions.
#'
#' @param regions1 A GRanges object
#' @param regions2 A GRanges object
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
