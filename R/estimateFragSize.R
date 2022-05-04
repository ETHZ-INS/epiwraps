#' estimateFragSize
#'
#' @param bam The path to one or more bam files
#' @param ctrl Optional path to a control bam file (if `length(bam)>1`, the same
#'   control will be used for all).
#' @param binSize Bin size. The precision of the reported fragment size is 
#'   necessary lower than this. A higher bin size will improve the summit 
#'   identification in low-coverage regions. We recommend leaving the default
#'   value.
#' @param mfold The range of fold-enrichment over the control (if `ctrl` 
#'   provided) or of coverages for the identification of regions based on which
#'   distance will be estimated.
#' @param minSummitCount The minimum read count for a summit to be considered.
#' @param useSeqLevels An optional vector of seqLevels in which to conduct the
#'   analysis.
#' @param maxSize The maximum size of regions to be used
#' @param priorLength The prior fragment length (use for read extension to 
#'   identify enriched regions)
#' @param ret The type of return, either a 'table' of pairs of summits and their
#'   properties, or a 'plot', or the median/mean/mode of the distance distribution.
#' @param blacklist Optional `GRanges` of blacklisted regions to be excluded.
#' @param BPPARAM A `BiocParallel` parameter object for multithreading. Only 
#'   used if multiple files are given in `bam`.
#'
#' @return By default, the estimated (mode) fragment length(s), but see the 
#'   `ret` argument
#' @export
#'
#' @examples
#' # get an example bam file
#' bam <- system.file("extdata", "ex1.bam", package="Rsamtools")
#' # create bigwig
#' estimateFragSize(bam)
estimateFragSize <- function(bam, ctrl=NULL, binSize=10L, mfold=c(10,50), ...,
                             minSummitCount=8L, useSeqLevels=NULL,
                             maxSize=2500L, priorLength=200L, blacklist=NULL, 
                             ret=c("mode", "median", "mean", "tables", "plots", "distances"),
                             BPPARAM=SerialParam()){
  ret <- match.arg(ret)

  if(is.list(bam)){
    if(!all(c("cov","cop","con","msize") %in% names(bam)))
      stop("Unrecognized input")
    # internal use by other functions, to avoid reading coverages again
    co <- bam
  }else if(length(bam)>1){
    # multiple bam files against a common input
    # we want to read the input only once
    stopifnot(all(unlist(lapply(bam,file.exists(bam)))))
    if(!is.null(ctrl)){
      if(!is(ctrl, "RleList")){
        if(length(ctrl)>1)
          stop("Only one `ctrl` bam can be provided. If you have multiple ",
               "pairs of signal and control bam files, please make independent",
               "calls to the function.")
        ctrl <- Reduce("+", bamChrChunkApply(ctrl, keepSeqLvls=useSeqLevels, 
                                             paired=FALSE, FUN=function(x){
                   x <- suppressWarnings(resize(x, pmax(width(x),priorLength)))
                   coverage(trim(x))
                 }))
      } 
    }
    ret2 <- ifelse(ret=="plots", "tables", ret)
    names(bam) <- bam
    res <- bplapply(bam, ctrl=ctrl, binSize=binSize, mfold=mfold, ret=ret2, 
                    maxSize=maxSize, useSeqLevels=useSeqLevels,
                    minSummitCount=minSummitCount, ...,
                    FUN=estimateFragSize, BPPARAM=BPPARAM)
    if(ret2!="tables") return(unlist(res))
    res <- lapply(setNames(names(res[[1]]),names(res[[1]])), FUN=function(x){
      dplyr::bind_rows(lapply(res, FUN=function(y) y[[x]]), .id="sample")
    })
    if(ret=="tables") return(res)
    return(.fragLengthPlots(res))
  }else{
    # obtain the coverages
    co <- .getStrandedCoverages(bam, keepSeqLvls=useSeqLevels, binSize=5L,
                                fragLength=priorLength)
  }
  if(!is.null(ctrl)){
    # convert coverage to foldchange
    if(!is(ctrl, "RleList")){
      ctrl <- Reduce("+", bamChrChunkApply(ctrl, keepSeqLvls=useSeqLevels, 
                                           paired=FALSE, FUN=function(x){
        x <- suppressWarnings(resize(x, pmax(width(x),priorLength)))
        coverage(trim(x))
      }))
    }
    nf <- .covTrimmedMean(co$cov)/.covTrimmedMean(ctrl)
    co$cov <- setNames((1L+co$cov)/(1L+ctrl*nf), names(co$cov))
  }
  msize <- co$msize
  # identify mfold-based peaks for estimation
  regions <- slice(co$cov, mfold[1], mfold[2], rangesOnly=TRUE)
  # keep only regions larger than the median read size
  regions <- GRanges(rep(factor(names(regions),names(regions)), lengths(regions)),
                     unlist(regions))
  regions <- regions[which(width(regions)>msize)]
  regions <- reduce(regions, min.gapwidth=msize)
  regions <- regions[width(regions)<=maxSize]
  
  # enlarge around the summit
  regions <- trim(resize(.viewl2gr(viewRangeMaxs(Views(co$cov, regions))), 
                         maxSize, fix="center"))
  # remove potential blacklist & overlaps
  if(!is.null(blacklist)) regions <- regions[overlapsAny(regions, blacklist)]
  
  if(ret %in% c("tables","plots"))
    d1 <- .strandedCovTable(regions, co)
  
  # for each region, find the + and - summits and their distance
  regions <- split(ranges(regions), seqnames(regions), drop=TRUE)
  names(seqlvls) <- seqlvls <- names(regions)
  dn <- dplyr::bind_rows(lapply(seqlvls, FUN=function(x){
    regions <- regions[[x]]
    sp <- .getSummits(co$cop[[x]], mfold, regions)
    sn <- .getSummits(co$con[[x]], mfold, regions)
    data.frame(peak.size=unlist(width(regions)),
               distance=end(sn)-start(sp),
               score.pos=as.numeric(mcols(sp)$score),
               score.neg=as.numeric(mcols(sn)$score))
  }))
  # difference between summit counts shouldn't be too large
  dn$absDiff <- abs(dn$score.pos-dn$score.neg)
  # keep only pairs of summits that are convincing
  dn <- dn[abs(dn$distance)<maxSize & abs(dn$distance)>msize & 
             dn$absDiff<(3*median(dn$absDiff)) &
             (dn$score.pos+dn$score.neg)>=minSummitCount ,]
  if(nrow(dn)<50)
    warning("A low number of regions was retained to estimate fragment size.",
            "You may consider increase the `mfold` range.")
  if(ret=="tables") return(list(perPeakDistance=dn, coverages=d1))
  if(ret=="mean") return(mean(abs(dn$distance)))
  if(ret=="median") return(median(abs(dn$distance)))
  if(ret=="distances") return(abs(dn$distance))
  if(ret=="mode"){
    d <- density(abs(dn$distance))
    return(round(d$x[which.max(d$y)]))
  }
  .fragLengthPlots(list(perPeakDistance=dn, coverages=d1))
}

.strandedCovTable <- function(regions, co){
  regions <- regions[start(regions)>=1]
  vp <- Views(co$cop, regions)
  vp <- Reduce("+",as(vp[lengths(vp)>0], "IntegerList"))
  vn <- Views(co$con, regions)
  vn <- lapply(vn, FUN=function(x) Reduce("+",as(x, "IntegerList")))
  vn <- Reduce("+",as(vn[lengths(vn)>0], "IntegerList"))
  size <- length(vp)
  data.frame(rel_pos=seq(from=-floor(size/2), to=ceiling(size/2), length.out=length(vp)),
             count=c(vp,vn), strand=rep(c("+","-"),each=length(vp)))
}

.fragLengthPlots <- function(x, span=0.05){
  if(!requireNamespace("ggplot2", quietly=TRUE))
    stop("The 'ggplot2' package is required.")
  dn <- x$perPeakDistance
  d1 <- x$coverages
  q <- quantile(abs(dn$distance), c(0.05,0.5,0.95))
  p1 <- ggplot(dn, aes(abs(distance))) + geom_histogram(bins=50) + 
    geom_vline(xintercept=q, linetype=c("dashed","solid","dashed")) +
    annotate("label", x=q[2], label=round(q[2]), y=2) +
    labs(x="Per-peak fragment length distribution", y="Count")
  p2 <- ggplot(d1, aes(rel_pos, count, colour=strand)) + geom_line(alpha=0.5) + 
    geom_smooth(formula=y~x, method="loess", span=span, size=1.2) + 
    labs(x="Relative position", y="Read start count")
  if(requireNamespace("cowplot", quietly=TRUE))
    return(plot_grid(p1, p2, nrow=2))
  print(p1)
  print(p2)
}

# sets coverages below a threshold to 0
.covListMin <- function(cov, minC=0L){
  if(minC<=0) return(cov)
  as(lapply(cov, FUN=function(x){
    runValue(x)[runValue(x)<minC] <- 0L
    x
  }), "RleList")
}

.getSummits <- function(cop, mfold=NULL, regions=NULL){
  if(is.null(mfold) & is.null(regions))
    stop("One of 'mfold' or 'regions' must be given.")
  if(is.null(regions)){
    sp <- slice(cop, lower=mfold[1], upper=mfold[2])
    g <- resize(viewRangeMaxs(sp), width=1L, fix="center")
  }else{
    sp <- Views(cop, regions)
    g <- resize(viewRangeMaxs(sp), width=1L, fix="center")
  }
  if(is(cop, "RleList"))
    g <- GRanges(rep(factor(names(g),names(g)), lengths(g)), unlist(g))
  mcols(g)$score <- max(sp)
  g
}

.getStrandedCoverages <- function(bam, keepSeqLvls=NULL, fragLength=100L, binSize=5L){
  covs <- bamChrChunkApply(bam, keepSeqLvls=keepSeqLvls, FUN=function(x){
    msize <- median(width(x))
    x <- trim(suppressWarnings(resize(x, width=pmax(width(x),fragLength))))
    cov <- coverage(x)
    x <- trim(resize(x,binSize))
    pos <- coverage(x[strand(x)=="+"])
    neg <- coverage(x[strand(x)=="-"])
    list(cov=cov, pos=pos, neg=neg, nreads=length(x), msize=msize)
  })
  list(
    msize=median(unlist(lapply(covs, FUN=function(x) x$msize)), na.rm=TRUE),
    cov=Reduce("+", lapply(covs, FUN=function(x) x$cov)),
    cop=Reduce("+", lapply(covs, FUN=function(x) x$pos)),
    con=Reduce("+", lapply(covs, FUN=function(x) x$neg))
  )
}