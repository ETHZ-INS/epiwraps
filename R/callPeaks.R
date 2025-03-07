#' callPeaks
#' 
#' This is a R-based implementation of the general MACS2 strategy (Zhang et al.,
#'  Genome Biology 2008), taking some freedom here and there in comparison to 
#'  the original. The function is still under development, especially with 
#'  respect to single-end reads, where some optimization might still be needed.
#'  For paired-end reads, the results are nearly identical with those of MACS2, 
#'  with two main differences: 1) the p-values are more conservative (and 
#'  arguably better calibrated) and 2) because the implementation does not rely 
#'  on sliding windows, with default settings the peaks are narrower.
#'
#' @param bam A signal bam file
#' @param ctrl An optional (but highly recommended) control bam file
#' @param paired Logical, whether the reads are paired
#' @param binSize Binsize used to estimate peak shift
#' @param fragLength Fragment length. Ignored if `paired=TRUE`. This is only 
#'   used for the initial candidate region identification, and sizes are 
#'   adjusted after, so it doesn't need to be very precise.
#' @param minPeakCount The minimum summit count for a region to be considered.
#'   Decreasing this can substantially increase the running time.
#' @param minFoldChange The minimum fold-change for a region to be considered.
#'   Decreasing this can substantially increase the running time.
#' @param blacklist An optional `GRanges` of regions to be excluded (or the path
#'   to such a file). Since the blacklisted regions are removed from both the
#'   signal and control peaks, this also has an important impact on the 
#'   empirical FDR (when `ctrl` is given).
#' @param bgWindow The windows to consider (in addition to the peak itself) 
#'   for local background.
#' @param type The type of peaks to identify ('narrow' or 'broad').
#' @param outFormat The output format ('custom' or 'narrowPeak')
#' @param pthres The p-value threshold to use
#' @param verbose Logical; whether to output progress messages
#' @param ... Passed to \code{\link{bamChrChunkApply}}
#'
#' @return A `GRanges`
#' 
#' @details 
#' `callPeaks` takes about twice as long to run as MACS2, and uses more memory.
#' If dealing with very large files (or a very low memory system), consider
#' increasing the number of processing chunks, for instance with `nChunks=10`.
#' 
#' The function uses \code{\link{bamChrChunkApply}} to obtain the coverages,
#' and can accept any argument of that function. This means that the 
#' `mapqFilter` and bam `flgs` arguments can be used to restrict the reads used.
#' 
#' @importFrom IRanges Views viewMaxs viewMeans slice viewRangeMaxs relist IntegerList
#' @importFrom stats setNames pnorm ppois optim density
#' @importFrom S4Vectors mean.Rle
#' @export
callPeaks <- function(bam, ctrl=NULL, paired=FALSE, type=c("narrow","broad"), 
                      nullH=c("local","global.nb","global.bin"), ### TO IMPLEMENT
                      blacklist=NULL, binSize=10L, fragLength=200L, 
                      minPeakCount=5L, minFoldChange=1.3, pthres=10^-3,
                      maxSize=NULL, bgWindow=c(1,5,10)*1000, pseudoCount=1L,
                      useStrand=TRUE, outFormat=c("custom", "narrowPeak"),
                      verbose=TRUE, ...){
  type <- match.arg(type)
  if(is.null(maxSize)) maxSize <- ifelse(type=="narrow",1000L,5000L)
  outFormat <- match.arg(outFormat)
  binSize <- as.integer(binSize)
  minPeakCount <- as.integer(minPeakCount)
  fragLength <- as.integer(fragLength)
  if(!is.null(blacklist) && is.character(blacklist)){
    if(verbose) message("Importing blacklist...")
    blacklist <- rtracklayer::import(blacklist)
  }
  
  if(is(bam, "RleList")){
    useStrand <- FALSE
    if(verbose) message("Identifying candidate regions...")
    o <- .cpGetCandidates(bam, paired=paired, ctrl=ctrl,  verbose=verbose,
                          binSize=binSize, fragLength=fragLength, bgWindow=bgWindow, 
                          minFoldChange=minFoldChange, minPeakCount=minPeakCount,
                          pseudoCount=pseudoCount, useStrand=useStrand,
                          maxSize=maxSize)
    o <- list(o)
  }else{
    if(verbose) message("Reading signal and identifying candidate regions...")
    o <- bamChrChunkApply(bam, ..., paired=paired, ctrl=ctrl, verbose=verbose,
                          binSize=binSize, fragLength=fragLength, bgWindow=bgWindow, 
                          minFoldChange=minFoldChange, minPeakCount=minPeakCount,
                          pseudoCount=pseudoCount, useStrand=useStrand,
                          breakPeaks=type=="narrow", FUN=.cpGetCandidates)
  }
  if(!is.null(ctrl)){
    # adjust the per-block normalization factors
    nf <- unlist(lapply(o, FUN=function(x) x$nf), use.names=FALSE)
    mnf <- median(nf)
    if((max((nf-mnf))/mnf)>0.3)
      warning("There are large variations in the per-block normalization ",
              "factor estimates. This can happen when the signal is heavily ",
              "biased towards some chromosomes, and can lead to some peaks ",
              "being lost. To be on the safe side, you could reduce the number",
              " of blocks")
    for(i in seq_along(nf)){
      o[[i]]$regions$bg <- o[[i]]$regions$bg*nf[[i]]/mnf
      o[[i]]$negR$bg <- o[[i]]$negR$bg*mnf/nf[[i]]
    }
    covtm <- median(unlist(lapply(o, FUN=function(x) x$covtm), use.names=FALSE))
    negR <- unlist(GRangesList(lapply(o, FUN=function(x) x$negR)), use.names=FALSE)
    if(!is.null(blacklist)) negR <- negR[!overlapsAny(negR, blacklist)]
    
  }
  r <- unlist(GRangesList(lapply(o, FUN=function(x) x$r)), use.names=FALSE)
  if(!is.null(blacklist)) r <- r[!overlapsAny(r, blacklist)]
  if(verbose) message("Identified ", length(r), " candidate regions")
  rm(o)
  gc(full=TRUE)
  
  if(verbose) message("Computing significance...")
  if(!is.null(ctrl)){
    # calculate new logFC
    r$log10FE <- as.integer(round(100*log10((pseudoCount+r$meanCount)/(pseudoCount+r$bg))))
    r <- r[which(r$log10FE > as.integer(floor(100*log10(minFoldChange))))]
    # poisson deviation from background
    r$log10p <- -log10(ppois(r$meanCount, r$bg, lower.tail=FALSE))
    
    if(verbose) message("Computing false discovery rate using negative peaks...")
    # call peaks in the control
    p <- -log10(ppois(negR$meanCount, pmax(covtm, pseudoCount, negR$bg),
                      lower.tail=FALSE))
    metadata(r)$ctrl.log10p <- p
    mcols(r) <- cbind(mcols(r), getEmpiricalFDR(r$log10p, p))
  }else{
    if(verbose)
      message("(In the absence of a control, FDR is unlikely to be calibrated)")
    if(type=="narrow"){
      p <- pmax(ppois(r$meanCount, pmax(r$bg,pseudoCount), lower.tail=FALSE),
                ppois(r$maxPos, pseudoCount, lower.tail=FALSE),
                ppois(r$maxNeg, pseudoCount, lower.tail=FALSE))
      fc <- (pseudoCount+r$meanCount)/(pseudoCount+r$bg)
    }else{
      p <- ppois(r$meanCount, pseudoCount, lower.tail=FALSE)
      fc <- (pseudoCount+r$meanCount)/pseudoCount
    }
    r$log10FE <- as.integer(round(100*log10(fc)))
    r$log10p <- -log10(p)
    r$log10FDR <- round(-log10(p.adjust(p, method="holm")),2)
  }
  r <- r[r$log10p >= -log10(pthres)]
  if(verbose) message("Reporting ", length(r), " regions, ",
                      sum(r$log10FDR>-log10(0.05))," with FDR<0.05")
  r$log10p <- round(r$log10p,2)
  
  r$score <- as.integer(round(1000*(r$log10FE/quantile(r$log10FE, .99))))
  r$score[r$score>1000L] <- 1000L
  if(outFormat=="NarrowPeak") r <- .customPeaks2NarrowPeaks(r)
  r
}



.cpGetCandidates <- function(x, ctrl=NULL, paired=FALSE, blacklist=NULL, 
                             binSize=5L, fragLength=300L, minPeakCount=5L, minSize=10L,
                             maxSize=2000L, minFoldChange=1.3, bgWindow=c(1,5,10)*1000, 
                             pseudoCount=1L, useStrand=TRUE,
                             flgs=scanBamFlag(), breakPeaks=TRUE, verbose=TRUE, ...){
  covtm <- nf <- negR <- fsq <- NULL
  if(!is(x, "RleList")){
    if(verbose) message("Reading signal coverage...")
    if(paired){
      fsq <- quantile(width(x), c(0,0.05,0.1,0.25,0.5,0.75,0.9,0.95,1))
    }else{
      x <- trim(suppressWarnings(resize(x, width=fragLength, fix="start")))
    }
    co <- coverage(x)
    if(useStrand){
      cop <- coverage(resize(x[which(as.factor(strand(x))=="+")], binSize, fix="start"))
      con <- coverage(resize(x[which(as.factor(strand(x))=="-")], binSize, fix="start"))
    }
    rm(x)
    gc(full=TRUE, verbose=FALSE)
  }
  r <- slice(co, lower=minPeakCount)
  r <- r[width(r)>=minSize]
  if(verbose) message(sum(lengths(r)), " initial candidate regions")
  v <- Reduce("c", lapply(r[which(lengths(r)>0L)], FUN=RleList))
  r <- .viewl2gr(r, maxCount=unlist(viewMaxs(r)),
                 meanCount=unlist(viewMeans(r)), cov=v)
  rm(v)
  if(!is.null(ctrl)){
    if(!is(ctrl, "RleList")){
      if(verbose) message("Reading control coverage...")
      lvls <- lengths(co)[lengths(runValue(co))>1L]
      p <- ScanBamParam(flag=flgs,
                        which=GRanges(names(lvls), IRanges(1L, width=lvls)))
      if(paired){
        ctrl <- coverage(GRanges(readGAlignmentPairs(ctrl, param=p)))
      }else{
        ctrl <- GRanges(readGAlignments(ctrl, param=p))
        ctrl <- suppressWarnings(resize(ctrl, width=fragLength, fix="start"))
        ctrl <- coverage(trim(ctrl))
      }
    }
    covtm <- .covTrimmedMean(ctrl)
    nf <- covtm/.covTrimmedMean(co)
    fc <- nf*S4Vectors::mean.Rle(r$cov)/unlist(viewMeans(Views(ctrl, r)))
    r <- r[which(fc>minFoldChange)]
  }
  if(breakPeaks){
    r2 <- r[width(r)>maxSize]
    r2 <- suppressWarnings({
      .breakRegions2(r2, v=r2$cov, minW=round(fragLength/2), maxW=maxSize)
    })
    v <- Views(co, r2)
    r2$maxCount <- unlist(viewMaxs(v))
    r2$meanCount <- unlist(viewMeans(v))
    r2$cov <- Reduce("c", lapply(v, FUN=RleList))
    r <- sort(c(r[which(width(r)<=maxSize)], r2))
  }
  if(!is.null(blacklist)) r <- r[!overlapsAny(r,blacklist)]
  if(useStrand){
    if(verbose) message("Getting strand information")
    v <- Views(cop,r)
    r$wPos <- unlist(start(resize(viewRangeMaxs(v), width=1L, fix="center")))
    r$pos <- Reduce("c", lapply(v, FUN=RleList))
    v <- Views(con,r)
    r$wNeg <- unlist(start(resize(viewRangeMaxs(v), width=1L, fix="center")))
    r$neg <- Reduce("c", lapply(v, FUN=RleList))
    r$maxPos <- max(r$pos)
    r$maxNeg <- max(r$neg)
    r <- r[which(r$maxPos >= floor(minPeakCount/2) & 
                   r$maxNeg >= floor(minPeakCount/2))]
    rm(cop, con)
    if(breakPeaks){
      r <- .refinePeaks(r)
      v <- Views(co, r)
      r$maxCount <- unlist(viewMaxs(v))
      r$meanCount <- unlist(viewMeans(v))
    }
    rm(v)
    r$wNeg <- r$wPos <- r$pos <- r$neg <- NULL
  }
  r$cov <- NULL
  if(length(r)==0){
    r$bg <- vector(mode = "integer", length = 0L)
  }else if(!is.null(ctrl)){
    if(verbose) message("Computing neighborhood background")
    r$bg <- .getLocalBackground(ctrl, gr=r, windows=bgWindow)/nf
    r$log10FE <- as.integer(round(100*log10((pseudoCount+r$meanCount)/(pseudoCount+r$bg))))
    r <- r[which(r$log10FE > as.integer(floor(100*log10(minFoldChange))))]
    if(verbose) message("Obtaining negative peaks")
    # get negative regions
    negR <- slice(ctrl, lower=max(minPeakCount, minPeakCount/nf))
    negR <- .viewl2gr(negR, maxCount=unlist(viewMaxs(negR)),
                      meanCount=unlist(viewMeans(negR)))
    negR$bg <- nf*.getLocalBackground(co, gr=negR, windows=bgWindow)
    negR$log10FE <- as.integer(round(100*log10((pseudoCount+negR$meanCount)/(pseudoCount+negR$bg))))
    negR <- negR[which(negR$log10FE > as.integer(floor(100*log10(minFoldChange))))]
  }else{
    r$bg <- .getLocalBackground(co, r, windows=max(bgWindow), incRegion=FALSE)
  }
  if(verbose) message(length(r), " peaks after filtering")
  rm(co,ctrl)
  gc(verbose=FALSE)
  list(regions=r, fsq=fsq, negR=negR, nf=nf, covtm=covtm)
}

# refines peaks based on stranded read starts
.refinePeaks <- function(r, f=20, minC=3, minW=NULL, maxW=NULL){
  if(is.null(r$wPos)){
    r$wPos <- which.max()
  }
  if(is.null(maxW) || is.null(minW)){
    x <- (r$wNeg-r$wPos)[which(r$maxPos>(2*minC) & r$maxNeg>(2*minC))]
    x <- x[x>0]
    if(is.null(maxW)) maxW <- as.integer(round(quantile(x, 0.95)))
    if(is.null(minW)) minW <- as.integer(round(max(0.9*quantile(x, 0.05),20)))
  }
  medpo <- median(r$pos, na.rm=TRUE)
  medneg <- median(r$neg, na.rm=TRUE)
  w0 <- r$maxPos > minC & r$maxPos > f*medpo & !is.infinite(r$wPos) &
         r$maxNeg > minC & r$maxNeg > f*medneg & !is.infinite(r$wNeg) &
         (r$wNeg-r$wPos) >= minW
  w <- which(w0)
  start(r)[w] <- r$wPos[w]
  end(r)[w] <- r$wNeg[w]
  if(length(w <- which(!w0 & width(r)>maxW))>0){
    wPos <- round(.rleMedWhich(r$pos[w]>medpo[w]))
    wNeg <- round(.rleMedWhich(r$neg[w]>medneg[w]))
    w2 <- which((wNeg-wPos)>=minW)
    start(r)[w][w2] <- start(r)[w][w2]+wPos[w2]-1L
    end(r)[w][w2] <- start(r)[w][w2]+wNeg[w2]-1L
  }
  trim(suppressWarnings(resize(r, pmax(minW,width(r)), fix="center")))
}

.rleMedWhich <- function(rle){
  median(cumsum(runLength(rle))[IRanges::which(runValue(rle))])
}

# very slow implementation, not used -- see .breakRegions2 instead
.breakRegions <- function(x, v, maxW=250L, minW=50L){
  x$V <- v
  x$iC <- max(v)
  toFrag <- width(x)>maxW
  notFrag <- x[which(!toFrag)]
  x <- x[toFrag]
  i <- 0L
  while(length(x)>0){
    w <- which(x$V==x$iC)
    l1 <- unlist(w)
    l2 <- rep(width(x), lengths(w))-l1
    w2 <- relist(l1>=minW & l2>=minW, w)
    bp <- w[relist(l1>=minW & l2>=minW, w)]
    fragged <- lengths(bp)>0
    if(sum(fragged)>0){
      bp <- bp[fragged]
      bp <- unlist(bp[setNames(splitAsList(ceiling(lengths(bp)/2), seq_along(bp)),
                        NULL)])
      fr1 <- fr2 <- x[fragged]
      fr1 <- trim(suppressWarnings(resize(x[fragged], bp-1L)))
      fr1$V <- fr1$V[IntegerList(lapply(bp,seq_len))]
      fr2 <- trim(suppressWarnings(resize(x[fragged], width(x[fragged])-bp-1L, fix="end")))
      fr2$V <- fr2$V[IntegerList(mapply(from=bp,to=lengths(fr2$V),FUN=seq,SIMPLIFY=FALSE))]
      notFrag <- c(notFrag, fr1[which(width(fr1)<=maxW)], fr2[which(width(fr2)<=maxW)])
      x <- c(x[which(!fragged)], fr1[which(width(fr1)>maxW)], fr2[which(width(fr2)>maxW)])
    }
    i <- i + 1L
    rv <- runValue(x$V)
    x$iC <- max(rv[rv<x$iC])
    if(any(wStop <- x$iC<1L)){
      notFrag <- c(notFrag, x[which(wStop)])
      x <- x[which(!wStop)]
    }
  }
  sort(notFrag)
}

# x is a gr, v is a coverage RleList corresponding to elements of x
.breakRegions2 <- function(x, v, minW=25L, maxW=1000L, denom=2){
  w <- width(x)>maxW
  xO <- granges(x[which(!w)])
  x <- x[which(w)]
  ll <- mapply(x=v, start=start(x), SIMPLIFY=FALSE, FUN=function(x, start){
    x <- slice(x,lower=max(x)/denom, rangesOnly=TRUE)
    shift(x, shift=start-1L)
  })
  ll <- IRangesList(ll)
  gr <- GRanges(rep(seqnames(x), lengths(ll)),
                ranges=unlist(ll, recursive=FALSE, use.names=FALSE))
  gr <- reduce(gr, min.gapwidth=round(minW/2))
  gr <- reduce(resize(gr, pmax(minW, width(gr)), fix="center"))
  seqinfo(gr) <- seqinfo(xO)
  gr <- trim(gr)
  sort(c(xO, gr))
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
  if(is.null(th)) th <- quantile(unlist(runValue(x)), q)
  if(is(x,"RleList"))
    return(weighted.mean(unlist(lapply(x, th=th, FUN=.covTrimmedMean)),
                         lengths(x), na.rm=TRUE))
  w <- which(runValue(x)<=th & runValue(x)>0) 
  x <- Rle(runValue(x)[w], lengths=runLength(x)[w])
  mean(x)
}


.getLocalBackground <- function(co, gr, windows=c(1,5)*1000, incRegion=TRUE){
  v <- Views(co, trim(suppressWarnings(resize(gr, max(windows), "center"))))
  m <- do.call(cbind, lapply(sort(windows, decreasing=TRUE), FUN=function(w){
    if(w!=unlist(width(v))[1])
      v <- resize(v, w, fix="center")
    unlist(viewMeans(v))
  }))
  if(incRegion) m <- cbind(unlist(viewMaxs(Views(co, gr))),m)
  rowMaxs(m)
}

# resize peaks using read starts
.resizePeaks <- function(gr, fragLength, within=NULL, pos=NULL, neg=NULL, minC=4L){
  fragLength <- quantile(abs(fragLength), c(0.05,0.95))
  if(!is.null(pos) && !is.null(neg)){
    gr2 <- resize(gr, pmax(width(gr), fragLength[1]), fix="center")
    v <- Views(pos, gr2)
    posSummit <- resize(viewRangeMaxs(v), 1L, fix="center")
    posSummitCount<- viewMaxs(v)
    v <- Views(neg, gr2)
    negSummit <- resize(viewRangeMaxs(v), 1L, fix="center")
    negSummitCount<- viewMaxs(v)
  }else{
    il <- IntegerList(r$pos)
    posSummit <- round(median(which(il==max(il))))
    posSummitCount <- max(il)
    il <- IntegerList(r$neg)
    negSummit <- round(median(which(il==max(il))))
    negSummitCount <- max(il)
  }
  w <- which( posSummitCount>=minC & negSummitCount>=minC &
                abs(posSummitCount-negSummitCount)<
                rowMeans(cbind(posSummitCount,negSummitCount)) &
                abs(posSummit-negSummit) >= min(fragLength))
  if(sum(w)>0){
    # use stranded peaks
    ss <- cbind(posSummit,negSummit)[w,]
    if(!is.null(pos) && !is.null(neg)){
      end(gr)[w] <- rowMaxs(ss)
      start(gr)[w] <- rowMins(ss)
    }else{
      end(gr)[w] <- start(gr)[w]+rowMaxs(ss)-1L
      start(gr)[w] <- start(gr)[w]+rowMins(ss)-1L
    }
  }
  w <- setdiff(which(width(gr)>max(fragLength)), w)
  if(length(w)>0 && !is.null(within)){
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
          unlist(as(x,"IRangesList"), use.names=FALSE), ...)
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
#'   p-value, eventually extrapolating and preserving the ranking.
#'
#' @param log10p -log10 p-values of candidate peaks
#' @param pneg -log10 p-values of negative peaks
#' @param n The number of hypotheses (used when the empirical FDR is zero)
#'
#' @return A data.frame with the empirical FDR and a smoothed -log10(FDR)
#' @importFrom stats ecdf loess predict p.adjust
getEmpiricalFDR <- function(log10p, pneg, n=length(log10p)*10){
  log10p[is.infinite(log10p)] <- max(log10p[!is.infinite(log10p)])
  pneg[is.infinite(pneg)] <- max(pneg[!is.infinite(pneg)])
  holm <- -log10(p.adjust(10^-log10p, method="holm", n=n))
  holm[is.infinite(holm)] <- max(holm[!is.infinite(holm)])
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
  # ensure that the log10 FDR is not larger than Holm's
  logFDR <- pmin(logFDR, holm, na.rm=TRUE)
  data.frame(empiricalFDR=eFDR, log10FDR=round(logFDR,2))
}

.makeMonotonic <- function(x, by){
  d <- data.frame(id=seq_along(x), x=x, by=by)
  d <- d[order(d$by),]
  d$x <- cummax(d$x)
  d[order(d$id), "x"]
}

