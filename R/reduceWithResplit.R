#' Merge regions, re-splitting large merges using local overlap minima
#'
#' This is an alternative to something like 
#' \code{reduce(unlist(GRangesList(peaks)))} for merging overlapping regions"
#' it tries to break up large merged regions based on the profile of overlap
#' with the un-merged regions. We typically use this to merge for instance peaks
#' called on different samples.
#'
#' @param peaks A list of \code{\link[GenomicRanges]{GRanges-class}}, or a 
#'   \code{\link[GenomicRanges]{GRanges-class}} containing overlapping regions.
#' @param softMaxSize The (merged) peak size below which re-splitting will be 
#'   attempted
#' @param relTroughDepth The minimum depth of local minima, as a fraction of the
#'   maximum. E.g. with a maxima of 12 peaks, the default of 1/4 would require 
#'   the minima to be below or equal to 9.
#' @param minTroughDepth The absolute minimum depth of local minima, in number 
#'   of peaks below the maxima.
#' @param minTroughWidth The minimum width of the local minima.
#' @param minDistFromBoundary The minimum distance of the local minima from the
#'   peak border.
#' @param minPeakSize The minimum final peak size.
#' @param BPPARAM BiocParallel Param object for multithreading. If set, 
#'   chromosomes are split into threads.
#'
#' @return A reduced `GRanges` of non-overlapping peaks.
#' 
#' @details
#' This is an alternative to something like 
#' \code{reduce(unlist(GRangesList(peaks)))}, which stitches overlapping regions
#' together and can result in large regions that can be problematic for some
#' applications. The function tries to break those large regions into composing
#' by using the coverage by the original (un-merged) regions. See the example 
#' below for an illustration.
#' The procedure first reduces `peaks`, then identifies reduced regions whose 
#' width is above a certain threshold (`softMaxSize`). For those regions, a 
#' coverage by the original peaks is computed to identify local minima
#' ('troughs') in the coverage that could divide the region into sub-regions of 
#' desirable lengths. `relThroughDepth` determines the minimum depth of the 
#' trough (i.e. decrease) as a fraction of the maximum coverage in the region,
#' while `minTroughDepth` determines the absolute minimum depth.
#' Note that the algorithm iterates through regions one by one and as such is 
#' quite slow, hence multithreading is recommended for large sets of regions.
#' 
#' @export
#' @importFrom BiocParallel SerialParam bplapply
#' @import GenomicRanges
#' @examples
#' # consider the following example set of regions:
#' gr <- GRanges("1", IRanges(c(100,120,140,390,410,430,120),
#'                            width=rep(c(200,520),c(6,1))))
#' plotSignalTracks(list(regions=gr, "# overlapping regions"=coverage(gr),
#'                       reduced=reduce(gr)), region=reduce(gr))
#' # if we are interested in having smaller regions, clearly it would seem 
#' # sensible here to cut roughly in the middle, since we have two distinct 
#' # groups of regions that are only joined by a single region
#' 
#' (redGr <- reduceWithResplit(gr, softMaxSize=100))
#' plotSignalTracks(list("# overlapping regions"=coverage(gr),
#'                       reduced=reduce(gr), "reduced\n\\w resplit"=redGr),
#'                       region=reduce(gr))
reduceWithResplit <- function(peaks, softMaxSize=500L, relTroughDepth=1/3, 
                              minTroughDepth=2L, minTroughWidth=1L,
                              minDistFromBoundary=150L, minPeakSize=100L,
                              BPPARAM=BiocParallel::SerialParam()){
  if(is.list(peaks) || is(peaks, "GRangesList"))
    peaks <- sort(unlist(GRangesList(peaks)))
  p <- reduce(peaks, with.revmap=TRUE)
  p1 <- p[width(p)<=softMaxSize | lengths(p$revmap)==1]
  p1$revmap <- NULL
  p <- p[width(p)>softMaxSize & lengths(p$revmap)>1]
  v <- Views(coverage(peaks), p)
  gaps <- bplapply(v, BPPARAM=BPPARAM, FUN=function(v){
    unlist(IRangesList(lapply(seq_along(v), FUN=function(i){
      x <- v[[i]]
      newR <- IRanges(1L,length(x))
      prevgaps <- gaps <- IRanges()
      curgaps <- 1L
      while(length(curgaps)>0 && any(width(newR)>softMaxSize)){
        curgaps <- unlist(IRangesList(lapply(seq_along(newR), FUN=function(j){
          x <- x[start(newR[j]):end(newR[j])]
          cth <- max(0,min(max(x)-ceiling(max(x)*relTroughDepth), max(x)-minTroughDepth))
          
          x <- .doSplitIR(x, cth, minDistFromBoundary=minDistFromBoundary,
                          minTroughWidth=minTroughWidth, minPeakSize=minPeakSize)
          shift(x, start(newR[j])-1L)
        })))
        gaps <- reduce(sort(c(gaps, curgaps)))
        if(identical(gaps,prevgaps)) curgaps <- NULL
        prevgaps <- gaps
        newR <- setdiff(IRanges(1L,length(x)), gaps)
      }
      shift(gaps, start(v)[i]-1L)
    })))
  })
  gaps <- as(IRangesList(gaps), "GRanges")
  sort(c(p1,setdiff(p, gaps)))
}

.doSplitIR <- function(x, cth, minDistFromBoundary=100L, minTroughWidth=10L, minPeakSize=100L){
  vs <- slice(x, upper=cth)
  ir <- as(vs, "IRanges")
  mcols(ir)$min <- min(vs)
  w <- start(ir)>=minDistFromBoundary & end(ir)<=length(x)-minDistFromBoundary
  ir <- ir[width(ir)>=minTroughWidth & 
             (!w | (length(x)-width(ir))>=minPeakSize)]
  w <- start(ir)>=minDistFromBoundary & end(ir)<=length(x)-minDistFromBoundary
  if(length(ir)==0) return(ir)
  if(length(which(w))>0){
    # there is a non-boundary window, use that
    # if many, use the deepest trough
    ir <- ir[mcols(ir)$min==min(mcols(ir)$min)]
    # if many, use the most central one
    w <- which.min(abs(length(x)/2-(start(ir)+width(ir)/2)))
    ir2 <- ir[w]
  }else{
    # boundary window
    ir <- ir[which.max(width(ir))]
    ir <- resize(ir, floor(width(ir)/2),
                 fix=ifelse(start(ir)>1,"start","end"))
    x2 <- x[start(ir):end(ir)]
    vs <- slice(x2, upper=min(x2))
    ir2 <- shift(as(vs, "IRanges"), start(ir)-1L)
    mcols(ir2)$min <- min(vs)
    w <- start(ir2)>=minDistFromBoundary & end(ir2)<=(length(x)-minDistFromBoundary)
    ir2 <- ir2[w]
    if(length(ir2)==0) return(ir2)
    ir2 <- ir2[which(mcols(ir2)$min==min(mcols(ir2)$min))]
    ir2 <- ir2[which.max(width(ir2))]
  }
  ir2 <- resize(ir2,pmin(width(ir2),minTroughWidth),fix="center")
  # last double-check
  ir2 <- ir2[start(ir2)>=minDistFromBoundary & end(ir2)<=(length(x)-minDistFromBoundary)]
  while(any(w <- (c(start(ir2),length(x))-c(0,end(ir2)[-length(ir2)]))<minPeakSize)){
    ir2 <- ir2[-which(w)[1]]
  }
  ir2
}
