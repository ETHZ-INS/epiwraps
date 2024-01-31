#' meltSignals
#'
#' Aggregates and melts a list of signal matrices, for plotting (with ggplot).
#'
#' @param ml A named list of matrices as produced by \code{\link{signal2Matrix}}
#' @param fun An optional custom aggregation function (or named list thereof).
#' @param trim The quantile above which to trim values. If a numeric vector of 
#'   length 2, will be used as lower and upper quantiles beyond which to trim.
#'
#' @return A data.frame.
#' @export
#' @importFrom matrixStats colMedians
meltSignals <- function(ml, fun=NULL, splitBy=NULL, trim=0.98){
  stopifnot(is.list(ml))
  ml <- .comparableMatrices(ml, checkAttributes=TRUE)
  ml <- .applyTrimming(ml, trim)
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
#' @param trim Quantiles trimmed at each extreme (for linear normalization)
#'
#' @return A renormalized list of signal matrices.
#' @export
#' @importFrom edgeR calcNormFactors
renormalizeBorders <- function(ml, method="linear", trim=NULL,
                               nWindows=max(floor(ncol(ml[[1]])/10),1)){
  ml <- .comparableMatrices(ml, checkAttributes=TRUE)
  b <- do.call(cbind, lapply(ml, FUN=function(x){
    as.numeric(cbind(x[,seq_len(nWindows)],
                     x[,seq(from=ncol(x)-nWindows+1, to=ncol(x))]))
  }))
  if(is.null(trim)) trim <- (sum(b!=0)/length(b))/10
  if(method=="linear"){
    nf <- apply(b, 2, FUN=function(x){
      y <- mean(x, trim=trim)
      if(y==0) y <- mean(x)
      y
    })
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

.resizeNmatrix <- function(x, i = NULL, j = NULL){
  a <- attributes(x)
  xcl <- class(x)
  a$dimnames[[1]] <- a$dimnames[[1]][i]
  a$dimnames[[2]] <- a$dimnames[[2]][j]
  a$names <- a$names[i]
  if(!is.null(i)) x <- x[i, , drop = FALSE]
  if(!is.null(j)) x <- x[, j, drop = FALSE]
  a$dim <- dim(x)
  attributes(x) <- a
  class(x) <- xcl
  x
}


#' clusterSignalMatrices
#'
#' @param ml A list of signal matrices, as produced by \code{\link{signal2Matrix}}.
#' @param k The number of clusters to generate
#' @param scaleRows Logical; whether to scale rows for clustering
#' @param scaleCols Logical; whether to scale columns (i.e. signals/samples)
#' @param assay Assay to use (ignored unless `ml` is an ESE object)
#' @param use What values to use for clustering. By default, uses 
#'   \code{\link[EnrichedHeatmap]{enriched_score}}. Other options are 'full' 
#'   (uses the full signal for clustering), 'max' (uses the maximum value in 
#'   the region), or 'center' (use the value at the center of the region).
#' @param by Optional factor/character/integer vector of the same length as 
#'   `ml`. When scaling rows, this can be used to indicate which rows should be
#'   scaled together.
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
                                  by=rep(1L,length(ml)),
                                  assay=1L, trim=c(0.05,0.95), nstart=3, ...){
  if(is(ml, "SummarizedExperiment")) ml <- .ese2ml(ml, assay=assay)
  #ml <- .comparableMatrices(ml)
  k <- unique(as.integer(k))
  stopifnot(all(k>1 & k<nrow(ml[[1]])))
  use <- match.arg(use)
  
  .getMatCenter <- function(x){
    a <- attributes(x)
    if(length(ti <- a$target_index)==0)
      ti <- c(max(a$upstream_index),min(a$downstream_index))
    rowMeans(x[,ti,drop=FALSE])
  }
  
  ml <- lapply(ml, FUN=function(x){
    q <- quantile(x, trim)
    x[x>q[2]] <- q[2]
    x[x<q[1]] <- q[1]
    if(scaleCols) x <- (x-mean(x))/sd(x)
    switch(use,
           full=x,
           max=matrixStats::rowMaxs(x),
           center=.getMatCenter(x),
           enrich=enriched_score(x))
  })
  if(scaleRows){
    stopifnot(length(ml)==length(by))
    ml <- lapply(split(ml, by), FUN=function(x){
      m <- do.call(cbind, x)
      m <- m - rowMeans(m)
      sd <- sqrt(matrixStats::rowVars(m, na.rm=TRUE))
      sd[which(sd==0)] <- 1
      m <- m/sd
    })
  }
  m <- do.call(cbind, ml)
  res <- lapply(setNames(k,k), FUN=function(x){
    cl <- kmeans(dist(m), centers=x)
    ve <- round(100*sum(cl$betweenss)/sum(c(cl$withinss,cl$betweenss)))
    list(cl=factor(as.character(cl$cluster),as.character(seq_len(x))), ve=ve)
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