#' meltSignals
#'
#' Aggregates and melts a list of signal matrices, for plotting (with ggplot).
#'
#' @param ml A named list of matrices as produced by \code{\link{signal2Matrix}}
#' @param fun An optional custom aggregation function (or named list thereof).
#'
#' @return A data.frame.
#' @export
#' @importFrom matrixStats colMedians
meltSignals <- function(ml, fun=NULL, splitBy=NULL){
  stopifnot(is.list(ml))
  ml <- .comparableMatrices(ml, checkAttributes=TRUE)
  if(!is.null(splitBy)){
    stopifnot(length(splitBy)==nrow(ml[[1]]))
    return(dplyr::bind_rows(lapply(split(seq_along(splitBy), splitBy),
                                   FUN=function(i){
      meltSignals(lapply(ml, i=i, FUN=.resizeNmatrix), fun=fun)
    }), .id="split"))
  }
  a <- attributes(ml[[1]])
  x <- c( -1*(a$extend[1]/length(a$upstream_index))*
            rev(seq_along(a$upstream_index)),
          (a$extend[2]/length(a$downstream_index))*
            seq_along(a$downstream_index) )
  d <- data.frame( row.names=NULL, position=rep(x, length(ml)),
                   sample=rep(factor(names(ml), names(ml)), each=length(x)) )
  d$mean <- unlist(lapply(ml, colMeans))
  d$SD <- unlist(lapply(ml, matrixStats::colSds))
  d$SE <- d$SD/sqrt(nrow(ml[[1]]))
  d$median <- unlist(lapply(ml, matrixStats::colMedians))
  if(!is.null(fun)){
    if(!is.list(fun)) fun <- list(fn.value=fun)
    stopifnot(all(unlist(lapply(fun,is.function))))
    if(is.null(names(fun)))
      names(fun) <- make.unique(rep("fn.value",length(fun)))
    for(fn in names(fun))
      d[[fn]] <- unlist(lapply(ml, FUN=function(x) apply(x,2,fun[[fn]])))
  }
  d
}

#' mergeSignalMatrices
#'
#' @param ml A named list of matrices as produced by \code{\link{signal2Matrix}}
#' @param aggregation The method to aggregate matrices
#'
#' @return A `normalizedMatrix` object
#' @export
mergeSignalMatrices <- function(ml, aggregation=c("mean","sum","median")){
  aggregation <- match.arg(aggregation)
  ml <- .comparableMatrices(ml, checkAttributes=TRUE)
  switch(aggregation,
    mean=Reduce("+", ml)/length(ml),
    sum=Reduce("+", ml),
    median=.medianSignalMat(ml)
  )
}

#' @importFrom matrixStats rowMedians
.medianSignalMat <- function(ml){
  x <- matrixStats::rowMedians(do.call(cbind,lapply(ml, as.numeric)))
  y <- ml[[1]]
  y[seq_len(nrow(y)),seq_len(ncol(y))] <- x
  y
}

#' renormalizeBorders
#'
#' This function renormalizes a list of signal matrices on the assumption that
#' the left/right borders of the matrices represent background signal which 
#' should be equal across samples.
#' \strong{This is not a safe normalization procedure}: it will work only if 
#' 1) the left/right borders of the matrices are sufficiently far from the 
#' signal (e.g. peaks), and 2) the signal-to-noise ratio is comparable across 
#' samples.
#'
#' @param ml A named matrix list as produced by \code{\link{signal2Matrix}}.
#' @param method Either "linear", or a normalization method, passed to 
#'   \code{\link[edgeR]{calcNormFactors}}.
#'
#' @return A renormalized list of signal matrices.
#' @export
#' @importFrom edgeR calcNormFactors
renormalizeBorders <- function(ml, method="linear", 
                               nWindows=max(floor(ncol(ml[[1]])/20),1)){
  ml <- .comparableMatrices(ml, checkAttributes=TRUE)
  b <- do.call(cbind, lapply(ml, FUN=function(x){
    as.numeric(cbind(x[,seq_len(nWindows)],
                     x[,seq(from=ncol(x)-nWindows+1, to=ncol(x))]))
  }))
  if(method=="linear"){
    nf <- apply(b, 2, trim=0.01, FUN=mean)
    nf <- nf/median(nf)
  }else{
    nf <- calcNormFactors(b, method=method, lib.size=rep(1,ncol(b)))
  }
  ml <- rescaleSignalMatrices(ml, 1/nf)
  ml
}

#' rescaleSignalMatrices
#'
#' @param ml A named matrix list as produced by \code{\link{signal2Matrix}}.
#' @param scaleFactors A numeric vector of same length as `ml`, 
#'   indicating the scaling factors by which to multiply each matrix.
#'   Alternatively, a numeric matrix with a number of rows equal to the length 
#'   of `ml`, and two columns indicating the alpha and beta arguments of a 
#'   s3norm normalization.
#'
#' @return A renormalized list of signal matrices.
#' @export
rescaleSignalMatrices <- function(ml, scaleFactors, applyLinearly=NULL){
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

.resizeNmatrix <- function(x, i=seq_len(nrow(x)), j=seq_len(ncol(x))){
  a <- attributes(x)
  xcl <- class(x)
  a$dimnames[[1]] <- a$dimnames[[1]][i]
  a$dimnames[[2]] <- a$dimnames[[2]][j]
  a$names <- a$names[i]
  x <- x[i,j,drop=FALSE]
  a$dim <- dim(x)
  attributes(x) <- a
  class(x) <- xcl
  x
}


#' bwNormFactors
#' 
#' Estimates normalization factors for a set of bigwig files.
#'
#' @param x A vector of paths to bigwig files, or alternatively a list of 
#'   coverages in RleList format.
#' @param wsize The size of the random windows. If any of the bigwig files 
#'  records point events (e.g. insertions) at high resolution (e.g. nucleotide),
#'  use a lot (e.g. >10k) of small windows (e.g. `wsize=10`), as per default
#'  settings. Otherwise the process can be lightened by using fewer bigger 
#'  windows.
#' @param nwind The number of random windows
#' @param peaks A list of peaks (GRanges) for each experiment in `x`, or a
#'   vector of paths to such files
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
#' signal-to-noise ratio. The 'MAnorm' approach (Shao et al., Genome Biology 
#' 2012) assumes that regions that are commonly enriched in two experiments
#' should on average have the same signal in the two experiments. These methods
#' then use linear scaling.
#' The 'S3norm' (Xiang et al., NAR 2020) and '2cLinear' methods try to 
#' normalize both simultaneously. S3norm does this in a log-linear fashion (as 
#' in the publication), while '2cLinear' does it on the original scale.
#'
#' @return
#' @export
bwNormFactors <- function(x, wsize=10L, nwind=20000L, peaks=NULL, trim=0.05,
                          useSeqLevels=NULL, 
                          method=c("background","SES","MAnorm","S3norm",
                                   "2cLinear")){
  method <- match.arg(method)
  chrsizes <- lapply(x, FUN=function(x){
    if(is.character(x)) return(seqlengths(BigWigFile(x)))
    if(is(x,"RleList")) return(lengths(x))
    stop("Unrecognized input")
  })
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

  if(method!="MAnorm"){
    windows <- .randomTiles(chrsizes, nwind, wsize)
    wc <- do.call(cbind, lapply(x, windows=windows, FUN=.getCovVals))
    wc <- wc[which(rowSums(is.na(wc))==0 & rowSums(wc)>0),]
    if(any(colSums(wc)<50))
      warning("Some samples have less than 50 non-zero windows. Consider ",
              "increasing the window size or (better) the number of windows.")
    if(method %in% c("background","SES")){
      nf <- apply(wc, 2, trim=trim, na.rm=TRUE, FUN=mean)
      return(setNames(median(nf, na.rm=TRUE)/nf, names(x)))
    }
  }

  stopifnot(length(peaks)==length(x))
  peaks <- lapply(peaks, FUN=function(x){
    if(is.character(x)) x <- rtracklayer::import(x)
    x
  })
  wc <- wc[!overlapsAny(windows, 
                        reduce(unlist(GRangesList(peaks), use.names=FALSE))),]
  
  cost <- function(p){
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
  
  nf <- lapply(seq_along(x)[-1], FUN=function(i){
    common <- reduce(resize(intersect(peaks[[1]],peaks[[i]]), width=wsize))
    common <- common[seqnames(common) %in% names(chrsizes)]
    common <- keepSeqlevels(common, seqlevelsInUse(common), pruning.mode="coarse")
    seqlengths(common) <- chrsizes[seqlevels(common)]
    x1 <- .getCovVals(x[[1]], common)
    x2 <- .getCovVals(x[[i]], common)
    if(method=="MAnorm") return(mean(x1,trim=trim)/mean(x2,trim=trim))
    optim(c(a=ifelse(method=="S3norm",1,0), b=1), fn=cost)$par
  })
  if(method=="MAnorm") return(setNames(c(1,unlist(nf)), names(x)))
  
  nf <- rbind(c(ifelse(method=="S3norm",1,0),1), do.call(rbind, nf))
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

.getCovVals <- function(x, windows){
  if(is.character(x)){
    x <- rtracklayer::import(x, format="BigWig", 
                             selection=BigWigSelection(windows))
    x <- coverage(x, weight=x$score)
  }
  windows <- sort(windows)
  y <- rep(NA_integer_, length(windows))
  w <- which(seqnames(windows) %in% names(x))
  if(length(w)==0) stop("No window found in coverage! Wrong seqlevel style?")
  windows <- windows[w]
  windows <- keepSeqlevels(windows, seqlevelsInUse(windows), pruning.mode="coarse")
  windows <- split(ranges(windows),seqnames(windows),drop=TRUE)
  y[w] <- unlist(viewMaxs(Views(x[names(windows)], windows)),
                 use.names=FALSE)
  y
}

#' clusterSignalMatrices
#'
#' @param ml A list of signal matrices, as produced by \code{\link{signal2Matrix}}.
#' @param k The number of clusters to generate
#' @param scaleRows Logical; whether to scale rows for clustering
#' @param scaleCols Logical; whether to scale columns (i.e. signals/samples)
#' @param use What values to use for clustering. By default, uses 
#'   \code{\link[EnrichedHeatmap]{enriched_score}}. Other options are 'full' 
#'   (uses the full signal for clustering), 'max' (uses the maximum value in 
#'   the region), or 'center' (use the value at the center of the region).
#' @param trim Values to trim (applied individidually for each signal matrix)
#' @param nstart Number of starts for k-means clustering
#' @param ... Passed to `kmeans`
#'
#' @return If `k` is of length 1, a vector of cluster labels, corresponding to 
#'   the rows of `ml`. If `length(k)>1`, a list of two data.frames containing:
#'   1) the cluster labels at the different resolutions, and 2) the variance 
#'   explained by clusters at each resolution.
#' @export
#' @importFrom matrixStats rowMaxs rowVars
#' @importFrom stats kmeans
clusterSignalMatrices <- function(ml, k, scaleRows=FALSE, scaleCols=FALSE,
                                  use=c("enrich","full","max","center"),
                                  trim=c(0.05,0.95), nstart=3, ...){
  ml <- .comparableMatrices(ml)
  k <- unique(as.integer(k))
  stopifnot(all(k>1 & k<nrow(ml[[1]])))
  use <- match.arg(use)
  ml <- lapply(ml, FUN=function(x){
    q <- quantile(x, trim)
    x[x>q[2]] <- q[2]
    x[x<q[1]] <- q[1]
    x
  })
  if(scaleCols) ml <- lapply(ml, FUN=function(x) (x-mean(x))/sd(x))
  ml <- switch(use,
               full=ml,
               max=lapply(ml, FUN=matrixStats::rowMaxs),
               center=lapply(ml, FUN=function(x){
                 a <- attributes(x)
                 if(length(ti <- a$target_index)==0)
                  ti <- c(max(a$upstream_index),min(a$downstream_index))
                 rowMeans(x[,ti,drop=FALSE])
               }),
               enrich=lapply(ml, enriched_score))
  m <- do.call(cbind, ml)
  if(scaleRows){
    m <- m - rowMeans(m)
    m <- m/sqrt(matrixStats::rowVars(m, na.rm=TRUE))
  }
  res <- lapply(setNames(k,k), FUN=function(x){
    cl <- kmeans(dist(m), centers=k)
    ve <- round(100*sum(cl$betweenss)/sum(c(cl$withinss,cl$betweenss)))
    list(cl=factor(as.character(cl$cluster),as.character(seq_len(k))), ve=ve)
  })
  if(length(res)==1){
    message("  ~", res[[1]]$ve, "% of the variance explained by clusters")
    return(res[[1]]$cl)
  }
  cl <- do.call(cbind, lapply(res, FUN=function(x) as.data.frame(x$cl)))
  colnames(cl) <- k
  ve <- data.frame(k=k, varExplained=unlist(lapply(res, FUN=function(x) x$ve)))
  list(clusters=cl, varExplained=ve)
}


# .mlSmoothedMax <- function(x, k=3){
#   if(is.list(x)){
#     x <- lapply(x, k=k, FUN=.mlSmoothedMax)
#     return(matrixStats::rowMaxs(do.call(cbind, x)))
#   }
#   x <- asplit(unclass(x), 1)
#   max(runmean(as(lapply(x, as.numeric), "RleList"), k))
# }