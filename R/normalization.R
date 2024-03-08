#' getNormFactors : estimate normalization factors from genomic signal files
#' 
#' Estimates normalization factors for a set of samples (i.e. bam/bigwig files).
#'
#' @param x A vector of paths to bigwig files, to bam files, or alternatively a 
#'   list of coverages in RleList format. (Mixed formats are not supported)
#' @param wsize The size of the random windows. If any of the bigwig files 
#'  records point events (e.g. insertions) at high resolution (e.g. nucleotide),
#'  use a lot (e.g. >10k) of small windows (e.g. `wsize=10`), as per default
#'  settings. Otherwise the process can be lightened by using fewer bigger 
#'  windows.
#' @param nwind The number of random windows
#' @param peaks A list of peaks (GRanges) for each experiment in `x`, or a
#'   vector of paths to such files, or a single GRanges of peaks to use for 
#'   MAnorm
#' @param trim Amount of trimming when calculating means
#' @param method Normalization method (see details)
#' 
#' @return A vector of normalization factors, or for the 'S3norm' and '2cLinear'
#'   methods, a numeric matrix with a number of rows equal to the length 
#'   of `x`, and two columns indicating the alpha and beta terms.
#' 
#' @details 
#' The 'background' or 'SES' normalization method (they are synonyms here)
#' (Diaz et al., Stat Appl Gen et Mol Biol, 2012) assumes that the background
#' noise should on average be the same across experiments, an assumption that 
#' works well in practice when there are not very large differences in 
#' signal-to-noise ratio. The implementation uses the trimmed mean number of 
#' reads in random windows with non-zero counts.
#' The 'MAnorm' approach (Shao et al., Genome Biology 
#' 2012) assumes that regions that are commonly enriched (i.e. common peaks) in 
#' two experiments should on average have the same signal in the two 
#' experiments. These methods then use linear scaling.
#' The 'S3norm' (Xiang et al., NAR 2020) and '2cLinear' methods try to 
#' normalize both simultaneously. S3norm does this in a log-linear fashion (as 
#' in the publication), while '2cLinear' does it on the original scale.
#'
#' @return A vector of normalization factors or, for methods 'S3norm' and 
#'   '2cLinear', a matrix of per-sample normalization parameters.
#' @export
#' @examples
#' # we get an example bigwig file, and use it twice:
#' bw <- system.file("extdata/example_atac.bw", package="epiwraps")
#' bw <- c(sample1=bw, sample2=bw)
#' # we indicate to use only chr1, because it's the only one in the example file
#' getNormFactors(bw, useSeqLevels="1")
getNormFactors <- function(x, wsize=10L, nwind=20000L, peaks=NULL, trim=0.01,
                          useSeqLevels=NULL, paired=NULL, ..., verbose=TRUE,
                          method=c("background","SES","MAnorm","S3norm",
                                   "2cLinear")){
  method <- match.arg(method)
  if(is(x, "GRanges") || length(x)==1)
    stop("Can't normalize a single sample on its own!")
  
  if(all(vapply(x, FUN=is.character, FUN.VALUE=logical(1)))){
    # input are file paths
    inType <- .parseFiletypeFromName(x, requireUnique=TRUE)
    stopifnot(inType %in% c("bam","bw"))
    if(inType=="bam"){
      if(is.null(paired))
        message("`paired` unspecified, assuming the data to be single-end...")
      paired <- FALSE
      chrsizes <- lapply(x, FUN=function(x) seqlengths(BamFile(x)))
    }else{
      chrsizes <- lapply(x, FUN=function(x) seqlengths(BigWigFile(x)))
    }
  }else if(all(vapply(x, class2="GRanges", FUN=is, FUN.VALUE=logical(1)))){
    # input are granges
    stop("Not currently implemented with this type of input.")
  }else if(all(vapply(x, class2="RleList", FUN=is, FUN.VALUE=logical(1)))){
    # input are RleList
    chrsizes <- lapply(x, lengths)
  }else{
    stop("`x` is of unknown format or contains multiple formats.")
  }
  
  tt <- table(unlist(lapply(chrsizes, names)))
  tt <- names(tt)[tt==length(x)]
  if(length(tt)==0) stop("The files appear to have no chromosome in common.")
  chrsizes <- do.call(cbind, lapply(chrsizes, FUN=function(x) x[tt]))
  chrsizes <- matrixStats::rowMins(chrsizes)
  names(chrsizes) <- tt
  if(!is.null(useSeqLevels)){
    if(!all(useSeqLevels %in% names(chrsizes)))
      stop("Some of the requested seqlevels are not found in the data!")
    chrsizes <- chrsizes[useSeqLevels]
  }
  stopifnot(length(chrsizes)>0)
  names(seqlvls) <- seqlvls <- names(chrsizes)
  
  if(!(method %in% c("background","SES")) && is.null(peaks))
    stop("The selected normalization method requires peaks.")
  if(is(peaks,"GRanges")) peaks <- list(peaks)
  peaks <- lapply(peaks, FUN=function(x){
    if(is.character(x)) x <- rtracklayer::import(x)
    if(!is.null(useSeqLevels)) x <- x[seqnames(x) %in% useSeqLevels]
    x
  })
  
  if(method=="MAnorm"){
    if(verbose) message("Comparing coverage in peaks...")
    return(.MAnorm(x, peaks, trim=trim, useSeqLevels=useSeqLevels, 
                   paired=paired, ...))
  }
  
  if(verbose) message("Comparing coverage in random regions...")
  windows <- .randomTiles(chrsizes, nwind, wsize)
  wc <- do.call(cbind, lapply(x, windows=windows, ..., FUN=.getCovVals))
  w <- which(rowSums(is.na(wc))==0 & rowSums(wc)>0)
  windows <- windows[w]
  wc <- wc[w,]
  if(any(colSums(wc)<50))
    warning("Some samples have less than 50 non-zero windows. Consider ",
            "increasing the window size or (better) the number of windows.")
  
  if(method %in% c("background","SES")){
    nf <- apply(wc, 2, FUN=function(x){
      if(sum(x!=0)<10)
        warning("Too few non-zero windows, please increase `nwind`.")
      y <- mean(x, trim=trim, na.rm=TRUE)
      if(y==0) y <- mean(x, na.rm=TRUE)
      y
    })
    return(setNames(median(nf, na.rm=TRUE)/nf, names(x)))
  }
  if(verbose) message("Comparing coverage in peaks...")
  wc <- wc[!overlapsAny(windows, 
                        reduce(unlist(GRangesList(peaks), use.names=FALSE))),]
  .S3norm(x, peaks, bgWindows=wc, chrsizes, trim=trim, method=method, 
          wsize=wsize, ...)
}

#' @describeIn getNormFactors deprecated in favor of getNormFactors
#' @export
bwNormFactors <- function(x, ...){
  .Deprecated(old="bwNormFactors", new="getNormFactors",
              msg=paste("bwNormFactors is deprecated, please gradually switch",
                        "to `getNormFactors`."))
  getNormFactors(x, ...)
}

.MAnorm <- function(x, peaks, trim=0.05, useSeqLevels=NULL, ...){
  if(length(peaks)==1) peaks <- peaks[[1]]
  stopifnot(is(peaks, "GRanges") || length(peaks)==length(x))
  if(is(peaks, "GRanges")){
    refP <- sort(peaks)
    peaks <- lapply(seq_len(length(x)), FUN=function(x) peaks)
    ref <- 1
  }else{
    ref <- .getRefSampleFromPeaks(peaks)
    refP <- sort(peaks[[ref]])
  }
  refP$ID <- seq_along(refP)
  peaks <- lapply(peaks, FUN=function(x){
    if(identical(refP,x)) return(NULL)
    refP[overlapsAny(refP, x)]$ID
  })
  refC <- .getCovVals(x[[ref]], refP, ...)
  nf <- mapply(p=peaks, x=x, FUN=function(p,x){
    if(is.null(p)) return(1)
    co <- .getCovVals(x, refP[p], ...)
    if(any(refC<0))
      return(mean(refC[p],trim=trim)/mean(co, trim=trim))
    nf <- edgeR::calcNormFactors(cbind(refC[p], co))
    nf[1]/nf[2]
  })
  setNames(nf, names(x))
}

.S3norm <- function(x, peaks, bgWindows, chrsizes, trim, wsize=10L, 
                    method="S3norm", ...){
  if(any(colSums(bgWindows)<50))
    warning("Some samples have less than 50 non-zero windows. Consider ",
            "increasing the window size or (better) the number of windows.")
  
  cost <- function(p, x1, x2, bg1, bg2){
    a <- p["a"]
    b <- p["b"]
    if(method=="S3norm"){
      x1 <- a*(x1)^b
      bg1 <- a*(bg1[bg1>0])^b
    }else{
      x1 <- a+b*x1
      bg1 <- a+(bg1[bg1>0])*b
    }
    abs(mean(x1, trim=trim)-mean(x2, trim=trim)) +
      abs(mean(bg1, trim=trim) - mean(bg2[bg2>0L], trim=trim))
  }
  
  common <- NULL
  if(length(peaks)==1){
    common <- peaks[[1]]
    ref <- 1L
  }else{
    ref <- .getRefSampleFromPeaks(peaks)
  }
  nf <- lapply(seq_along(x), FUN=function(i){
    if(is.null(common)){
      common <- reduce(resize(intersect(peaks[[ref]],peaks[[i]]), width=wsize))
      common <- common[seqnames(common) %in% names(chrsizes)]
      common <- keepSeqlevels(common, seqlevelsInUse(common), pruning.mode="coarse")
      seqlengths(common) <- chrsizes[seqlevels(common)]
    }
    x1 <- .getCovVals(x[[ref]], common, ...)
    x2 <- .getCovVals(x[[i]], common, ...)
    bg1 <- bgWindows[,ref]
    bg2 <- bgWindows[,i]
    optim(c(a=ifelse(method=="S3norm",1,0), b=1), 
          x1=x1, x2=x2, bg1=bg1, bg2=bg2, fn=cost)$par
  })
  
  nf <- do.call(rbind, nf)
  row.names(nf) <- names(x)
  colnames(nf) <- letters[1:2]
  attributes(nf)$applyLinearly <- method=="2cLinear"
  nf
}


#' @importFrom IRanges IRangesList
.randomTiles <- function(chrsizes, nwind, wsize){
  winPerChr <- round(nwind*(chrsizes/sum(chrsizes)))
  winPerChr <- winPerChr[winPerChr>=1]
  windows <- as(IRangesList(lapply(setNames(names(winPerChr),names(winPerChr)),
       FUN=function(x){
         possibleWindows <- floor(chrsizes[x]/wsize)
         IRanges(sort(1L+wsize*(sample.int(possibleWindows, winPerChr[x])-1L)), 
                 width=wsize)
       })), "GRanges")
  seqlengths(windows) <- chrsizes[seqlevels(windows)]
  windows
}

.getCovVals <- function(x, windows, paired=FALSE, ...){
  if(is.character(x)){
    if(.parseFiletypeFromName(x)=="bam"){
      if(paired){
        x <- readGAlignmentPairs(BamFile(x),
                                 param=ScanBamParam(which=windows, ...))
      }else{
        x <- readGAlignments(BamFile(x), param=ScanBamParam(which=windows, ...))
      }
      x <- coverage(x)
    }else{
      x <- rtracklayer::import(x, format="BigWig", as="RleList",
                               selection=BigWigSelection(windows))
    }
  }
  windows <- sort(windows)
  y <- rep(NA_integer_, length(windows))
  w <- which(as.factor(seqnames(windows)) %in% names(x))
  if(length(w)==0) stop("No window found in coverage! Wrong seqlevel style?")
  windows <- windows[w]
  windows <- keepSeqlevels(windows, seqlevelsInUse(windows), pruning.mode="coarse")
  windows <- split(ranges(windows),seqnames(windows),drop=TRUE)
  y[w] <- unlist(viewMaxs(Views(x[names(windows)], windows)),
                 use.names=FALSE)
  y
}

.getRefSampleFromPeaks <- function(peaks){
  if(is.list(peaks) && length(peaks)==1) return(1L)
  po <- sapply(seq_along(peaks), FUN=function(i){
    sapply(seq_along(peaks), FUN=function(j){
      if(i==j) return(NA_integer_)
      sum(overlapsAny(peaks[[i]], peaks[[j]]))
    })
  })
  ref <- which.max(matrixStats::rowMaxs(po,na.rm=TRUE))
  if(min(po[ref,],na.rm=TRUE)<100){
    if(min(po[ref,],na.rm=TRUE)<1){
      stop("Some of the experiments have no peak in common!")
    }
    warning("Some of the experiments have few (<100) peaks in common.",
            "The calculated factors might be inaccurate.")
  }
  ref
}

#' @export
#' @describeIn renormalizeSignalMatrices deprecated > renormalizeSignalMatrices
renormalizeBorders <- function(ml, trim=NULL, assay="input", nWindows=NULL){
  .Deprecated(
    old="renormalizeBorders", new="renormalizeSignalMatrices",
    msg=paste('renormalizeBorders is deprecated, please gradually',
              'switch to `renormalizeSignalMatrices(..., method="border"`.'))
  renormalizeSignalMatrices(ml, trim=trim, assay=assay, nWindows=nWindwos,
                            method="border")
}


#' renormalizeSignalMatrices
#'
#' Renormalizes a list of signal matrices or an EnrichmentSE object.
#'
#' @param ml A named matrix list or EnrichmentSE object as produced by 
#'   \code{\link{signal2Matrix}}.
#' @param method Either "border" or "top" (see details below).
#' @param trim Quantiles trimmed at each extreme before calculating 
#'   normalization factors.
#' @param fromAssay Assay to use (ignored unless `ml` is an EnrichmentSE 
#'   object), defaults to the first assay.
#' @param toAssay Assay in which to store the normalized data (ignored unless 
#'   `ml` is an EnrichmentSE object). By default an assay name will be set based
#'   on the normalization method used.
#' @param scaleFactors A numeric vector of same length as `ml`, 
#'   indicating the scaling factors by which to multiply each matrix.
#'   Alternatively, a numeric matrix with a number of rows equal to the length 
#'   of `ml`, and two columns indicating the alpha and beta arguments of a 
#'   s3norm normalization. Ignored unless `method="manual"`.
#'
#' @details
#' * `method="border"` works on the assumption that the left/right borders of the 
#' matrices represent background signal which should be equal across samples. As
#' a result, it will work only if 1) the left/right borders of the matrices are 
#' sufficiently far from the signal (e.g. peaks) to be chiefly noise, and 
#' 2) the signal-to-noise ratio is comparable across tracks/samples.
#' * `method="top"` instead works on the assumption that the highest signal should
#' be the same across tracks/samples.
#' By default, extreme values are trimmed before establishing either kind of 
#' normalization factor. The proportion trimmed can be set using the `trim` 
#' argument, and is by default 10% of the non-zero values.
#' * `method="manual"` enables the use of independently computed normalization
#' factors, for instance obtained through \code{\link{getNormFactors}}.
#'
#' @return Either a renormalized list of signal matrices or, if `ml` was an 
#'   `EnrichmentSE` object, the same object with an additional normalized
#'   assay automatically put at the front.
#' @export
#' @examples
#' # we first get an EnrichmentSE object:
#' bw <- system.file("extdata/example_atac.bw", package="epiwraps")
#' regions <- rtracklayer::import(system.file("extdata/example_peaks.bed", 
#'                                            package="epiwraps"))
#' m <- signal2Matrix(c(sample1=bw, sample2=bw), regions)
#' # we normalize them
#' m <- renormalizeSignalMatrices(m, method="border")
#' # see the `vignette("multiRegionPlot")` for more info on normalization.
renormalizeSignalMatrices <- function(ml, method=c("border","top","manual"), 
                                      trim=NULL, fromAssay="input", toAssay=NULL,
                                      nWindows=NULL, scaleFactors=NULL, ...){
  if(missing(method) && !is.null(scaleFactors)) method <- "manual"
  method <- match.arg(method)
  if(is(ml, "EnrichmentSE")){
    ml2 <- renormalizeSignalMatrices(.ese2ml(ml,assay=fromAssay), method=method,
                                     trim=trim, nWindows=nWindows,
                                     toAssay=toAssay, scaleFactors=scaleFactors)
    if(is.null(toAssay))
      toAssay <- switch(method,
                         "border"="borderNormalized",
                         "top"="topNormalized",
                         "normalized")
    return(addAssayToESE(ml, a=ml2, name=toAssay, replace=TRUE))
  }
  if(method=="manual"){
    stopifnot(!is.null(scaleFactors) && length(scaleFactors)==length(ml))
  }else{
    if(!is.null(scaleFactors)) 
      message('`scaleFactors` is ignored when method=!"manual"')
  }
  if(is.null(nWindows)) nWindows <- max(floor(ncol(ml[[1]])/10),1)
  ml <- .comparableMatrices(ml, checkAttributes=TRUE)
  b <- do.call(cbind, lapply(ml, FUN=function(x){
    x <- unclass(x)
    as.numeric(cbind(x[,seq_len(nWindows)],
                     x[,seq(from=ncol(x)-nWindows+1L,to=ncol(x))]))
  }))
  if(is.null(trim)) trim <- (sum(b!=0)/length(b))/10
  if(method=="linear"){
    nf <- apply(b, 2, FUN=function(x){
      y <- mean(x, trim=trim)
      if(y==0) y <- mean(x)
      y
    })
    nf <- nf/median(nf)
  }else if(method=="top"){
    nf <- apply(b, 2, FUN=function(x){
      y <- quantile(x, 1-trim)
      if(y==0) y <- max(x)
      y
    })
    nf <- nf/median(nf)
  }else if(method=="manual"){
    nf <- 1/scaleFactors
  }
  ml <- .applyMatricesScaling(ml, 1/nf)
  ml
}

# @param ml A named matrix list as produced by \code{\link{signal2Matrix}}.
# @param scaleFactors A numeric vector of same length as `ml`, 
#   indicating the scaling factors by which to multiply each matrix.
#   Alternatively, a numeric matrix with a number of rows equal to the length 
#   of `ml`, and two columns indicating the alpha and beta arguments of a 
#   s3norm normalization.
# @return A renormalized list of signal matrices.
.applyMatricesScaling <- function(ml, scaleFactors, applyLinearly=NULL){
  ml <- .comparableMatrices(ml, checkAttributes=TRUE)
  if(is.matrix(scaleFactors)){
    stopifnot(ncol(scaleFactors)==2 && nrow(scaleFactors)==length(ml))
    if(is.null(applyLinearly)){
      if(is.null(applyLinearly <- attributes(scaleFactors)[["applyLinearly"]]))
        stop("Could not determine whether the scaling factors should be",
             "applied linearly or on the log scale. Please use the ",
             "`applyLinearly` argument.")
    }
    for(i in seq_along(ml)){
      if(applyLinearly){
        ml[[i]] <- scaleFactors[i,1]+ml[[i]]*scaleFactors[i,2]
      }else{
        ml[[i]] <- scaleFactors[i,1]*ml[[i]]^scaleFactors[i,2]
      }
    }
  }else{
    stopifnot(is.numeric(scaleFactors) && length(scaleFactors)==length(ml))
    for(i in seq_along(scaleFactors)) ml[[i]] <- ml[[i]]*scaleFactors[i]
  }
  ml
}
