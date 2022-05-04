#' callPeaks
#' 
#' This is a R-based implementation of the general MACS2 strategy (Zhang et al.,
#'  Genome Biology 2008), taking some freedom here and there in comparison to 
#'  the original.
#'
#' @param bam A signal bam file
#' @param ctrl An optional (but highly recommended) control bam file
#' @param paired Logical, whether the reads are paired
#' @param binSize Binsize used to estimate peak shift
#' @param fragLength Fragment length. Ignored if `paired=TRUE`. If not given,
#'   will be estimated from the data.
#' @param minPeakCount The minimum summit count for a region to be considered.
#'   Decreasing this can substantially increase the running time.
#' @param minFoldChange The minimum fold-change for a region to be considered.
#'   Decreasing this can substantially increase the running time.
#' @param blacklist An optional `GRanges` of regions to be excluded (or the path
#'   to such a file). Since the blacklisted regions are removed from both the
#'   signal and control peaks, this also has an important impact on the 
#'   empirical FDR (when `ctrl` is given).
#' @param priorFragLength Prior fragment length (ignored if `paired=TRUE`). This
#'   is used for coverage computation to estimate peak shift.
#' @param bgWindow The windows to consider (in addition to the peak itself) 
#'   for local background.
#' @param type The type of peaks to identify ('narrow' or 'broad').
#' @param outFormat The output format ('custom' or 'narrowPeak')
#' @param pthres The p-value threshold to use
#' @param verbose Logical; whether to output progress messages
#' @param ... Passed to `estimateFragSize`
#'
#' @return A `GRanges`
#' 
#' @importFrom IRanges Views viewMaxs viewMeans slice viewRangeMaxs
#' @export
callPeaks <- function(bam, ctrl=NULL, paired=FALSE, type=c("narrow","broad"), 
                      blacklist=NULL, binSize=5L, fragLength=NULL, 
                      minPeakCount=5L, minFoldChange=1.3, pthres=10^-3, 
                      priorFragLength=150L, bgWindow=c(1,5,10)*1000, pseudoCount=1L,
                      outFormat=c("custom", "narrowPeak"), verbose=TRUE, ...){
  type <- match.arg(type)
  outFormat <- match.arg(outFormat)
  binSize <- as.integer(binSize)
  minPeakCount <- as.integer(minPeakCount)
  priorFragLength <- as.integer(priorFragLength)
  if(!is.null(fragLength)){
    fragLength <- as.integer(fragLength)
    stopifnot(fragLength>1)
    priorFragLength <- fragLength
  }
  if(!is.null(blacklist) && is.character(blacklist)){
    blacklist <- rtracklayer::import(blacklist)
  }
  
  if(!is.null(ctrl) && !is(ctrl, "RleList")){
    if(verbose) message("Reading control coverage...")
    # control (full) coverage
    if(paired){
      ctrl <- bamChrChunkApply(ctrl, paired=TRUE, FUN=coverage)
    }else{
      ctrl <- Reduce("+", bamChrChunkApply(ctrl, paired=FALSE, FUN=function(x){
        x <- suppressWarnings(resize(x, pmax(width(x), priorFragLength)))
        coverage(trim(x))
      }))
    }
  }
  
  # signal coverages
  if(verbose) message("Reading signal coverage...")
  if(paired){
    # we don't need to estimate fragment size, so just read full coverage
    co <- bamChrChunkApply(bam, paired=TRUE, FUN=function(x){
      list(fl=median(width(x)), cov=coverage(x))
    })
    fragLength <- median(unlist(lapply(co, FUN=function(x) x$fl)), na.rm=TRUE)
    co <- list(cov=Reduce("+", lapply(co, FUN=function(x) x$cov)))
  }else{
    # full + stranded coverages
    co <- .getStrandedCoverages(bam, binSize=binSize, fragLength=priorFragLength)
    if(is.null(fragLength)){
      if(verbose) message("Estimating fragment size...")
      fragLength <- estimateFragSize(co, ctrl, binSize=binSize, ...,
                                     priorFragLength=priorFragLength,
                                     blacklist=blacklist, ret="distances")
      if(verbose) message("  mean: ", round(mean(fragLength)),"; 90% CI: ",
                          paste(round(quantile(fragLength,c(0.05,0.95))),
                                collapse="-"))
    }
  }
  
  if(verbose) message("Identifying candidate regions...")
  bg <- .covTrimmedMean(co$cov)
  if(!is.null(ctrl)){
    # normalize ctrl
    nf <- bg/.covTrimmedMean(ctrl)
    ctrl <- setNames(ctrl*nf, names(ctrl))
    # calc foldchanges
    co$fc <- .rleFC(co$cov, ctrl)
    # identify regions of interest
    r1 <- slice(co$fc, lower=round(100*minFoldChange), rangesOnly=TRUE)
    r2 <- slice(co$cov, lower=minPeakCount, rangesOnly=TRUE)
    r <- Views(co$cov, intersect(r1, r2))
  }else{
    minPeakCount <- floor(max(minPeakCount,
                              .covTrimmedMean(co$cov)*minFoldChange))
    r <- slice(co$cov, lower=as.integer(minPeakCount))
  }
  
  if(type=="broad"){
    # merge nearby regions
    r <- reduce(r, min.gapwidth=round(median(fragLength)))
  }else{
    #r <- .breakPeaks(r, size=median(fragLength))
  }
  
  # GRanges of regions with stats
  rmax <- as(viewRangeMaxs(r) ,"GRanges")
  r <- .viewl2gr(r, summitCount=as.integer(unlist(viewMaxs(r))),
                 peakMean=round(unlist(viewMeans(r)),2),
                 summit=start(resize(rmax, 1L, fix="center")))
  if(!is.null(blacklist)) r <- r[!overlapsAny(r, blacklist)]

  # enlarge regions that are smaller than the 5% frag size
  q <- round(quantile(abs(fragLength),c(0.05,0.9)))
  r <- resize(r, pmax(width(r), q[1]), fix="center")
  if(verbose) message("  ", length(r), " candidate regions")
  
  if(type=="narrow" && !paired){
    # For regions that are too large, resize around the summit, keeping within
    # the regions' original boundaries
    if(length(w <- which(width(r)>q[2]))>0)
      ranges(r[w]) <- .resizeWithin(ranges(rmax[w]), q[2], within=r[w])
  }
  
  if(!is.null(ctrl)){
    message("Computing neighboring background and significance...")
    # calculate neighborhood background
    r$bgMax <- round(.getLocalBackground(ctrl, r, windows=bgWindow),2)
    # calculate new logFC
    r$logFC <- round(log2((pseudoCount+r$summitCount)/(pseudoCount+r$bgMax)),2)
    r <- r[r$logFC > log2(minFoldChange)]
    # poisson deviation from background
    r$log10p <- -log10(ppois(r$peakMean, r$bgMax, lower.tail=FALSE))
    message("Computing false discovery rate using the control...")
    # call peaks in the control
    fc2 <- .rleFC(ctrl, co$cov)
    fc2 <- slice(fc2, lower=as.integer(100*minFoldChange), rangesOnly=TRUE)
    if(sum(lengths(fc2))>1000)
      fc2 <- intersect(fc2, slice(ctrl, lower=minPeakCount, rangesOnly=TRUE))
    if(verbose) message("  ", length(fc2), " negative peaks")
    fc2 <- .viewl2gr(fc2)
    fc2 <- resize(fc2, pmax(width(fc2), mean(q)), fix="center")
    fc2$bg <- .getLocalBackground(co$cov, fc2, windows=bgWindow)
    fc2$cnt <- as.numeric(unlist(viewMeans(Views(ctrl, fc2)), use.names=FALSE))
    if(!is.null(blacklist)) fc2 <- fc2[!overlapsAny(fc2, blacklist)]
    p <- -log10(ppois(fc2$cnt, pmax(.covTrimmedMean(ctrl),pseudoCount,fc2$bg),
                      lower.tail=FALSE))
    metadata(r)$ctrl.log10p <- p
    mcols(r) <- cbind(mcols(r), getEmpiricalFDR(r$log10p, p))
  }else{
    message("Computing significance...")
    if(type=="narrow"){
      bg <- .getLocalBackground(co$cov, r, windows=max(bgWindow), incRegion=FALSE)
      p <- ppois(r$peakMean, bg, lower.tail=FALSE)
      r$logFC <- (pseudoCount+r$peakMean)/(pseudoCount+bg)
    }else{
      p <- ppois(r$peakMean, pseudoCount+.covTrimmedMean(co$cov), lower.tail=FALSE)
      r$logFC <- (pseudoCount+r$peakMean)/(pseudoCount+bg)
    }
    r$log10p <- -log10(p)
    r$log10FDR <- round(-log10(p.adjust(p, method="holm")),2)
  }
  r <- r[r$log10p >= -log10(pthres)]
  r$log10p <- round(r$log10p,2)
  r$score <- as.integer(round(1000*(r$logFC/quantile(r$logFC, .99))))
  r$score[r$score>1000L] <- 1000L
  if(outFormat=="NarrowPeak") r <- .customPeaks2NarrowPeaks(r)
  r
}

.customPeaks2NarrowPeak <- function(r){
  r$name <- Rle(rep(factor("."),length(r)))
  r$pointSource <- r$summit-start(r)
  mcols(r) <- mcols(r)[,c("name","score","logFC","log10p","log10FDR",
                          "pointSource")]
  r
}

# get trimmed non-zero mean of a Rle/RleList
.covTrimmedMean <- function(x, q=0.98, th=NULL){
  if(is.null(th)) th <- as.integer(quantile(unlist(runValue(x)), q))
  if(is(x,"RleList"))
    return(median(unlist(lapply(x, th=th, FUN=.covTrimmedMean), "RleList"),
                  na.rm=TRUE))
  w <- which(runValue(x)<=th & runValue(x)>0) 
  x <- Rle(runValue(x)[w], lengths=runLength(x)[w])
  mean(x)
}


.getLocalBackground <- function(co, gr, windows=c(1,5)*1000, incRegion=TRUE){
  v <- Views(co, resize(gr, max(windows), "center"))
  m <- do.call(cbind, lapply(sort(windows, decreasing=TRUE), FUN=function(w){
    if(w!=unlist(width(v))[1]) v <- resize(v, w, fix="center")
    unlist(viewMeans(v))
  }))
  if(incRegion) m <- cbind(unlist(viewMaxs(Views(co, gr))),m)
  rowMaxs(m)
}

# resize peaks using read starts
.resizePeaks <- function(gr, fragLength, within, pos=NULL, neg=NULL, minC=4L){
  fragLength <- quantile(abs(fragLength), c(0.05,0.95))
  gr2 <- resize(gr, pmax(width(gr), fragLength[1]), fix="center")
  v <- Views(pos, gr2)
  posSummit <- resize(viewRangeMaxs(v), 1L, fix="center")
  posSummitCount<- viewMaxs(v)
  v <- Views(neg, gr2)
  negSummit <- resize(viewRangeMaxs(v), 1L, fix="center")
  negSummitCount<- viewMaxs(v)
  w <- which( posSummitCount>=minC & negSummitCount>=minC &
                 abs(posSummitCount-negSummitCount)<
                   rowMeans(cbind(posSummitCount,negSummitCount)) &
                 abs(posSummit-negSummit) >= min(fragLength))
  if(sum(w)>0){
    # use stranded peaks
    ss <- cbind(posSummit,negSummit)[w,]
    end(gr)[w] <- rowMax(ss)
    start(gr)[w] <- rowMins(ss)
  }
  w <- setdiff(which(width(gr)>max(fragLength)), w)
  if(length(w)>0){
    # resize and try to keep within fc boundaries
    summits <- IRanges(gr$summit[w]-1L, gr$summit[w]+1L)
    ranges(gr)[w] <- .resizeWithin(summits, within=gr[w], size=min(fragLength))
  }
  gr
}


# resizes regions to a given size, but minimizing overflow from broader regions
.resizeWithin <-function(gr, size, within, doReduce=FALSE){
  size <- as.integer(size)
  if(doReduce){
    gr2 <- resize(gr, size, fix="center")
  }else{
    gr2 <- resize(gr, pmax(width(gr), size), fix="center")
  }
  wS <- start(gr2)<start(within)
  wE <- end(gr2)>end(within)
  if(length(w <- which(width(gr2)>width(within)))>0){
    # overflow on both sides
    ce <- rowMeans(cbind((start(gr2)[w]+start(within)[w])/2,
                         (end(gr2)[w]+end(within)[w])/2))
    end(gr2)[w] <- ceiling(ce+size/2)
    start(gr2)[w] <- floor(ce-size/2)
  }
  if(length(w <- which(width(gr2)<=width(within) & (wS | wE)))>0){
    # overflow on one side
    start(gr2)[w] <- pmax(start(gr2)[w], start(within)[w])
    end(gr2)[w] <- start(gr2)[w] + size
  }
  gr2
}

.viewl2gr <- function(x, ...){
  GRanges(rep(factor(names(x), names(x)), lengths(x)),
          unlist(as(x,"IRangesList")), ...)
}

.rleFC <- function(num, den, pseudoCount=1L, fact=100L, toInt=TRUE){
  if(is(num, "RleList"))
    return(as(mapply(num=num, den=den, pseudoCount=pseudoCount, toInt=toInt,
                     fact=fact, FUN=.rleFC), "RleList"))
  stopifnot(is(num, "Rle") && (is(den, "Rle") || is.numeric(den)))
  num <- fact*(pseudoCount+num)/(pseudoCount+den)
  if(toInt) runValue(num) <- as.integer(round(runValue(num)))
  num
}

#' getEmpiricalFDR
#' 
#' Computes the FDR as the proportion of negative peaks with a more extreme 
#'   p-value.
#'
#' @param log10p -log10 p-values of candidate peaks
#' @param pneg -log10 p-values of negative peaks
#' @param n The number of hypotheses (used when the empirical FDR is zero)
#'
#' @return A data.frame with the empirical FDR and a smoothed -log10(FDR)
#' @importFrom stats ecdf loess predict p.adjust
getEmpiricalFDR <- function(log10p, pneg, n=length(log10p)*10){
  holm <- -log10(p.adjust(10^-log10p, method="holm", n=n))
  # first, get the number >= extreme
  nInNeg <- (1-ecdf(pneg)(log10p))*length(pneg)
  nInPos <- (1-ecdf(log10p)(log10p))*length(log10p)
  eFDR <- nInNeg/(nInNeg+nInPos)
  logFDR <- -log10(eFDR)
  if(length(w <- which(is.na(logFDR) | is.infinite(logFDR)))>0){
    # replace FDRs of zero with Holm's, reducing the gap
    gap <- max(0,min(holm[w])-max(1+logFDR[-w]))
    logFDR[w] <- pmax(holm[w]-gap,2)
  }
  logFDR <- .makeMonotonic(logFDR, by=log10p)
  # smoothen
  logFDR <- predict(loess(logFDR~log10p, na.action="na.omit"), newdata=log10p)
  # ensure that the log10 FDR is not larger than Holm's
  logFDR <- pmin(logFDR, holm, na.rm=TRUE)
  data.frame(empiricalFDR=eFDR, logFDR=round(logFDR,2))
}

.makeMonotonic <- function(x, by){
  d <- data.frame(id=seq_along(x), x=x, by=by)
  d <- d[order(d$by),]
  d$x <- cummax(d$x)
  d[order(d$id),"x"]
}