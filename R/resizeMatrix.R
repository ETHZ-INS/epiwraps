#' resize a numeric matrix to given dimensions
#'
#' @param mat A numeric matrix
#' @param ndim The desired output dimensions
#' @param method Whether to use normal interpolation (`method="mean"`, 
#'   the default), or the max or min of the overlapping grid cells.
#'
#' @details 
#' For most cased this is based on Vyha's implementation (taken from 
#'  https://stackoverflow.com/a/23429527 ), but adding the possibility to 
#'  replace normal interpolation with the max/min of the overlapping grid cells.
#' This works well when the desired dimensions are at least half of the input
#' ones. When the desired dimensions are smaller, a different binning method 
#' is used, first applied on columns and then on rows.
#'
#' @return A numeric matrix of dimensions `ndim`
#' @importFrom stats approx
#' @export
resizeMatrix <- function(mat, ndim=dim(mat),
                         method=c("mean","max","min")){
  method <- match.arg(method)
  if(ndim[2]<(ncol(mat)/2)){
    mat <- Reduce(cbind, lapply(split(seq_len(ncol(mat)),
                                      cut(seq_len(ncol(mat)), ndim[2], 
                                          labels=FALSE)),
                           FUN=function(i){
                             x <- mat[,i,drop=FALSE]
                             switch(method,
                                    "mean"=rowMeans(x, na.rm=TRUE),
                                    "max"=matrixStats::rowMaxs(x, na.rm=TRUE),
                                    "min"=matrixStats::rowMins(x, na.rm=TRUE))
                           }))
  }
  if(ndim[1]<(nrow(mat)/2)){
    mat <- Reduce(rbind, lapply(split(seq_len(nrow(mat)),
                                      cut(seq_len(nrow(mat)), ndim[1], 
                                          labels=FALSE)),
                          FUN=function(i){
                             x <- mat[i,,drop=FALSE]
                             switch(method,
                                    "mean"=colMeans(x, na.rm=TRUE),
                                    "max"=matrixStats::colMaxs(x, na.rm=TRUE),
                                    "min"=matrixStats::colMins(x, na.rm=TRUE))
                           }))
  }
  if(all(dim(mat)==ndim)) return(mat)
  
  # input object
  obj <- list(x=seq_len(nrow(mat)), y=seq_len(ncol(mat)), z=mat)
  
  # output object
  ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])
  ndim <- dim(ans)
  
  # rescaling
  rescale <- function(x, newrange=range(x)){
    xrange <- range(x)
    mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
    newrange[1]+(x-xrange[1])*mfac
  }
  ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
  loc <- ncord
  loc[,1] <- rescale(ncord[,1], c(1,nrow(mat)))
  loc[,2] <- rescale(ncord[,2], c(1,ncol(mat)))
  
  # interpolation
  ans[ncord] <- .interp.surface(obj, loc, method=method)
  
  ans
}

# adapted from fields::interp.surface, adding alternatives to interpolation
.interp.surface <- function(obj, loc, method=c("mean","max","min")){
  method <- match.arg(method)
  x <- obj$x
  y <- obj$y
  z <- obj$z
  nx <- length(x)
  ny <- length(y)
  lx <- approx(x, 1:nx, loc[, 1])$y
  ly <- approx(y, 1:ny, loc[, 2])$y
  lx1 <- as.integer(floor(lx))
  ly1 <- as.integer(floor(ly))
  ex <- lx - lx1
  ey <- ly - ly1
  if(method %in% c("min","max")){
    fn <- pmax
    if(method=="min") fn <- pmin
    return(fn( z[cbind(lx1, ly1)], z[cbind(lx1 + as.integer(ex>0), ly1)],
               z[cbind(lx1, ly1 + as.integer(ey>0))],
               z[cbind(lx1 + as.integer(ex>0), ly1 + as.integer(ey>0))] ))
  }
  ex[lx1 == nx] <- 1
  ey[ly1 == ny] <- 1
  lx1[lx1 == nx] <- nx - 1
  ly1[ly1 == ny] <- ny - 1
  return(z[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) + 
           z[cbind(lx1 + 1, ly1)] * ex * (1 - ey) + 
           z[cbind(lx1, ly1 + 1)] * (1 - ex) * ey + 
           z[cbind(lx1 + 1, ly1 + 1)] * ex * ey)
}