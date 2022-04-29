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
estimateFragSize <- function(bam, ctrl=NULL, binSize=5L, mfold=c(10,50), ...,
                             minSummitCount=8L, useSeqLevels=NULL,
                             maxSize=2500L, priorLength=200L, blacklist=NULL, 
                             ret=c("mode", "median", "mean", "table", "plot"),
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
    ret2 <- ifelse(ret=="plot", "table", ret)
    names(bam) <- bam
    res <- bplapply(bam, ctrl=ctrl, binSize=binSize, mfold=mfold, ret=ret2, 
                    maxSize=maxSize, useSeqLevels=useSeqLevels,
                    minSummitCount=minSummitCount, ...,
                    FUN=estimateFragSize, BPPARAM=BPPARAM)
    if(ret2!="table") return(unlist(res))
    res <- dplyr::bind_rows(res, .id="sample")
    if(ret=="table") return(res)
    return(.plotDistDist(res))
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
    ctrl <- setNames(ctrl*sum(mean(co$cov))/sum(mean(ctrl)), names(ctrl))
    co$cov <- setNames((1L+co$cov)/(1L+ctrl), names(co$cov))
  }
  msize <- co$msize
  # identify mfold-based peaks for estimation
  regions <- slice(co$cov, mfold[1], mfold[2], rangesOnly=TRUE)
  # keep only regions larger than the median read size, and enlarge a bit
  regions <- GRanges(rep(factor(names(regions),names(regions)), lengths(regions)),
                     unlist(regions))
  regions <- regions[which(width(regions)>msize)]
  regions <- trim(resize(regions, width(regions)+as.integer(msize),
                         fix="center"))
  # remove potential blacklist & overlaps
  if(!is.null(blacklist)) regions <- regions[overlapsAny(regions, blacklist)]
  regions <- reduce(regions, min.gapwidth=msize)
  regions <- regions[width(regions)<=maxSize]
  # for each region, find the + and - summits and their distance
  regions <- split(ranges(regions), seqnames(regions), drop=TRUE)
  names(seqlvls) <- seqlvls <- names(regions)
  dn <- dplyr::bind_rows(lapply(seqlvls, FUN=function(x){
    regions <- regions[[x]]
    sp <- .getSummits(co$cop[[x]], mfold, regions)
    sn <- .getSummits(co$con[[x]], mfold, regions)
    data.frame(peak.size=unlist(width(regions)),
               distance=end(sp)-start(sn),
               score.pos=as.numeric(mcols(sp)$score),
               score.neg=as.numeric(mcols(sn)$score))
  }))
  # difference between summit counts shouldn't be too large
  dn$absDiff <- abs(dn$score.pos-dn$score.neg)
  # keep only pairs of summits that are convincing
  dn <- dn[abs(dn$distance)<maxSize & dn$absDiff<(3*median(dn$absDiff)) &
             (dn$score.pos+dn$score.neg)>=minSummitCount ,]
  if(nrow(dn)<50)
    warning("A low number of regions was retained to estimate fragment size.",
            "You may consider increase the `mfold` range.")
  if(ret=="table") return(dn)
  if(ret=="mean") return(mean(abs(dn$distance)))
  if(ret=="median") return(median(abs(dn$distance)))
  if(ret=="mode"){
    d <- density(abs(dn$distance))
    return(round(d$x[which.max(d$y)]))
  }
  .plotDistDist(dn)
}

# plots a histogram of the distances
.plotDistDist <- function(dn){
  dn <- dn[dn$distance<quantile(dn$distance, 0.98),]
  hist(abs(dn$distance), breaks=100, 
       main="Estimated fragment length", xlab="Length")
  abline(v=median(abs(dn$distance)), col="blue", lwd=3)
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