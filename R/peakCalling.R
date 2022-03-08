#' callPeaks
#' 
#' This is a very simple peak calling function wrapping around several 
#' \code{\link{csaw}} functions.
#'
#' @param signalBam Path to the ChIP bam file
#' @param controlBam Optional path to the input (i.e. negative control) bam file
#' @param broad Logical, whether to call broad regions rather than narrow peaks
#' @param fragSize Fragment length
#' @param wWidth Window width. Defaults to twice the fragment length 
#'   (`fragSize`) for narrow peaks and four times for broad.
#' @param spacing Shift between sliding windows
#' @param localSize Local neighborhood size
#' @param minLFC Minimum log2-foldchange
#' @param prior.count Prior count added to compute foldchange
#' @param filter The minimum count for a window to be considered
#' @param blacklist An optional GRanges of regions to be ignored.
#' @param p.threshold P-value threshold (poisson probability)
#' @param BPPARAM BiocParallel BPPARAM
#' @param verbose Logical; whether to print progress messages
#'
#' @return A GRanges object
#' 
#' @details 
#' The function supports broad and narrow peak calling, as well as calls based
#' on the comparison to a control (based on csaw's 
#' \code{\link[csaw]{filterWindowsControl}}), based on a local background for 
#' narrow peaks (see \code{\link[csaw]{filterWindowsLocal}}), or based on a 
#' global background for broad peaks (see 
#' \code{\link[csaw]{filterWindowsGlobal}}).
#'
#' @importFrom csaw windowCounts findMaxima regionCounts filterWindowsLocal
#' @importFrom csaw filterWindowsControl readParam
#' @importFrom SummarizedExperiment rowRanges rowRanges<- assay
#' @importFrom GenomicRanges width resize
#' @export
callPeaks <- function(signalBam, controlBam=NULL, broad=FALSE, fragSize=125L, 
                      wWidth=NULL, spacing=NULL, localSize=NULL, 
                      minLFC=log2(3), prior.count=3L, filter=1L+prior.count, 
                      p.threshold=10^-4, blacklist=NULL, BPPARAM=NULL, 
                      verbose=TRUE){
  stopifnot(is.logical(broad) && length(broad)==1)
  if(is.null(wWidth)) wWidth <- fragSize*2*(1+as.numeric(broad))
  if(is.null(spacing)) spacing <- max(20,round(wWidth/5))
  if(is.null(localSize)) localSize <- wWidth*4
  stopifnot(localSize>wWidth)
  stopifnot(prior.count>0)
  if(is.null(blacklist)) blacklist <- GRanges()
  rp <- readParam(minq=20, dedup=TRUE, discard=blacklist)
  if(is.null(controlBam)){
    if(verbose) message("Getting coverage across sliding windows...")
    wc <- windowCounts(signalBam, ext=fragSize, width=wWidth, filter=filter, 
                       param=rp, spacing=spacing)
    if(!broad){
      if(verbose) message("Evaluating neighborhoods...")
      wider <- suppressWarnings(resize(rowRanges(wc), localSize, fix="center"))
      wider <- regionCounts(signalBam, regions=wider, ext=fragSize, param=rp)
      filter.stat <- filterWindowsLocal(wc, wider, prior.count=prior.count)
    }else{
      if(verbose) message("Evaluating background...")
      bins <- windowCounts(signalBam, width=10000, bin=TRUE, param=rp)
      filter.stat <- filterWindowsGlobal(wc, bins, prior.count=prior.count)
    }
    rowRanges(wc)$count <- rowSums(assay(wc))
    rowRanges(wc)$background <- rowSums(assay(wider))*
      pmax(fragSize,width(wc))/width(wider)
    rm(wider)
    w <- which(filter.stat$filter>=minLFC)
    wc <- rowRanges(wc)[w]
    if(length(wc)==0) return(wc)
    score(wc) <- filter.stat$filter[w]
    wc$p.value <- ppois(wc$count, wc$background+prior.count, 
                        lower.tail=FALSE)
    wc <- wc[which(wc$p.value<p.threshold)]
  }else{
    if(verbose) message("Normalizing signal and input...")
    bins <- windowCounts(c(signalBam,controlBam), width=10000, bin=TRUE, 
                         param=rp, BPPARAM=BPPARAM)
    scale.info <- scaleControlFilter(bins[,1], bins[,2])
    rm(bins)
    if(verbose) message("Getting coverage across sliding windows...")
    wc <- windowCounts(c(signalBam,controlBam), ext=fragSize, width=wWidth, 
                       bin=TRUE, filter=filter, param=rp, BPPARAM=BPPARAM)
    if(verbose) message("Filtering...")
    filter.stat <- filterWindowsControl(wc[,1], wc[,2], prior.count=prior.count, 
                                        scale.info=scale.info)
    w <- which(filter.stat$filter>=minLFC)
    wc <- wc[w,1]
    rowRanges(wc)$count <- assay(wc)[,1]
    wc <- rowRanges(wc)
    if(length(wc)==0) return(wc)
    score(wc) <- filter.stat$filter[w]
    wc$p.value <- ppois(2^filter.stat$filter[w])    
  }
  wc <- wc[which(wc$p.value<p.threshold)]
  if(!broad){
    if(verbose) message("Removing overlapping peaks...")
    wc <- .removeOverlappingRanges(wc[order(-score(wc))])
  }
  re <- reduce(wc, with.revmap=TRUE, ignore.strand=TRUE)
  ag <- aggregate(as.data.frame(mcols(wc)[unlist(re$revmap),c("score","p.value")]), 
                  by=list(i=rep(seq_along(re), lengths(re$revmap))), FUN=max)
  mcols(wc) <- ag[,-1]
  wc
}


# Taken from scanMiR
.removeOverlappingRanges <- function(x, minDist=1L, ignore.strand=TRUE){
  red <- GenomicRanges::reduce(x, with.revmap=TRUE, min.gapwidth=minDist,
                               ignore.strand=ignore.strand)$revmap
  red <- red[lengths(red)>1]
  if(length(red)==0){
    if(retIndices) return(c())
    return(x)
  }
  i <- seq_along(x)
  toRemove <- c()
  while(length(red)>0){
    ## for each overlap set, we flag the index (relative to i) of the maximum
    ## (i.e. lowest in the list)
    top <- min(red) ## indexes of the top entry per overlap set, relative to i
    ## overlap of non-top entries to the top entries:
    o <- IRanges::overlapsAny(x[i[-top]],x[i[top]],maxgap=minDist)
    torem <- i[-top][which(o)] ## entries to remove, relative to x
    toRemove <- c(toRemove, torem) ## relative to x
    i <- setdiff(i,torem)
    ## and check again overlaps among this subset (revmap ind are relative to i)
    red <- GenomicRanges::reduce(x[i], with.revmap=TRUE, min.gapwidth=minDist,
                                 ignore.strand=ignore.strand)$revmap
    red <- red[lengths(red)>1]
  }
  if(length(toRemove)>0) x <- x[-toRemove]
  sort(x)
}