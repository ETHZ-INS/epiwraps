#' computeDeviationsFaster
#' 
#' A faster version of \code{\link[chromVAR]{computeDeviations}} (see details).
#'
#' @param counts A matrix of read counts per region(rows)/sample(columns), or a
#'   SummarizedExperiment with this as first assay, as produced by 
#'   \code{\link{getCounts}}. Can also be already normalized if 
#'   `normalize=FALSE`.
#' @param motifMatches A matrix of motif matches (either logical or with 
#'   something akin to binding probabilities), or a SummarizedExperiment 
#'   containing such an assay.
#' @param backgrounds A matrix of indices indicating background peaks, as 
#'   produced by \code{\link[chromVAR]{getBackgroundPeaks}}.
#' @param normalize Logical; whether to perform column sum normalization on 
#'  `counts`
#' @param welford Logical; whether to use Welford's online algorithm for 
#'   computing background means and standard deviations. If NULL (default), the 
#'   function will use Welford's if the predicted memory usage is above 10GB.
#' @param verbose Logical; whether to print messages, in particular the 
#'   projected memory size for large datasets.
#' @param BPPARAM An optional BiocParallel BPPARAM object for multi-threading.
#'
#' @return A SummarizedExperiment
#' 
#' @details
#' The results of this function are equivalent (within a very small precision) 
#' to those of \code{\link[chromVAR]{computeDeviations}}, but faster due to its
#' reliance on matrix operations, which enables the use of a larger number of
#' background iterations.
#' The use of matrix operations however implies that the entire results of each 
#' background iteration is stored in memory. While this is typically not a 
#' problem for bulk data, it often is with large (e.g. single-cell) datasets. 
#' For this reason, the function includes a variant using Welford's 
#' algorithm for the computation of the means and standard deviations. This is 
#' still faster than the original chromVAR implementation, but much less so than
#' when using `welford=FALSE`. By default, Welford's algorithm will be used if 
#' the projected memory usage is above 10GB.
#' The speed gain depends a lot on the dimensions of the different inputs and 
#' the number of threads used. The runtime when using Welford's algorithm can 
#' vary from 50 to 100% of the original chromVAR implementation. When not using 
#' Welford's algorithm, running times will typically be much smaller (typically 
#' roughly 20% of the original).
#' 
#' @author Pierre-Luc Germain
#' @export
computeDeviationsFaster <- function(counts, motifMatches, backgrounds,
                                    normalize=TRUE, welford=NULL, verbose=TRUE,
                                    BPPARAM=SerialParam(progress=TRUE)){
  stopifnot(nrow(counts)==nrow(motifMatches) &&
              nrow(counts)==nrow(backgrounds))
  stopifnot(min(b)>=1 & max(b)<=nrow(b))
  if(!is.integer(b[1,1]))
    stop("`backgrounds` should be a matrix of integers.")
  
  expectedMem <- .computeDeviationsMemoryUsage(counts, motifMatches, 
                                               backgrounds, BPPARAM$workers)
  if(is.null(welford)){
    welford <- (expectedMem[[1]]>expectedMem[[2]] && expectedMem[[1]]>=10)
    if(verbose && welford)
      message("Using Welford's online algorithm due to dataset size.\n",
              "This can be turned off with `welford=FALSE` to increase speed,",
              " but that may lead to very high memory usage.")
  }
  if(verbose){
    os <- ifelse(welford, expectedMem[[2]], expectedMem[[1]])
    if(os>0.2){
      os <- ifelse(os<0.1, paste0(round(os*1000,1), " MB"),
                   paste0(round(os,2), " GB"))
      message("Projected memory usage:", os)
    }
  }
  
  CD <- NULL
  if(inherits(counts, "SummarizedExperiment") ||
     inherits(counts, "SingleCellExperiment")){
    CD <- colData(counts)
    if("logcounts" %in% assayNames(counts)){
      if(verbose) message("Using the object's logcounts assay.")
      counts <- assay(counts, "logcounts")
      normalize <- FALSE
    }else{
      if(verbose && assayNames(counts)[1] != "counts"){
        message("Assuming the first assay of `counts` to be counts...")
      }
      counts <- assay(counts)
    }
  }
  if(normalize){
    # normalize the counts:
    counts <- t(10000*t(counts)/colSums(counts))
  }
  
  if(inherits(motifMatches, "SummarizedExperiment"))
    motifMatches <- assay(motifMatches)
  if(any(motifMatches>1))
    warning("motifMatches should be either binary or weights from 0 to 1.")
  
  # Ensure that the reordering is done on the least costly matrix
  doVariant <- ((is(mi, "sparseMatrix")-is(counts, "sparseMatrix"))+
                  (ncol(motifMatches)>ncol(counts))) > 0
  if(doVariant && is(counts, "CsparseMatrix")){
    counts <- as(counts, "RsparseMatrix")
  }else if(doVariant && is(motifMatches, "CsparseMatrix")){
    motifMatches <- as(motifMatches, "RsparseMatrix")
  }
  
  # computing background devs:
  if(welford){
    res <- bplapply(split(seq_len(ncol(backgrounds)), BPPARAM$workers),
                    BPPARAM=BPPARAM, FUN=function(idxs){
                      .computeBgDevChunkWelford(idxs, backgrounds, motifMatches,
                                                counts, doVariant=doVariant)
                    })
    # mean across bg iterations:
    m <- Reduce("+",lapply(res, \(x) x[[1]]))/length(res)
    # aggregate variances across bg sets and compute joint SDs
    sds <- sqrt(Reduce("+", lapply(res, \(x) x[[2]]^2)) + 
                  Reduce("+", lapply(res, \(x) (x[[1]]-m)^2)))
  }else{
    # computing background devs:
    res <- bplapply(seq_len(ncol(backgrounds)), BPPARAM=BPPARAM, 
                    FUN=function(i){
                      if(doVariant){
                        o <- t(crossprod(counts[backgrounds[,i],],motifMatches))
                      }else{
                        o <- crossprod(motifMatches[backgrounds[,i],],counts)
                      }
                      expect <- rowMeans(o)
                      (o-expect)/expect
                    })
    # mean and sd across bg iterations:
    m <- Reduce("+",res)/length(res)
    sds <- sqrt(Reduce("+",lapply(res, \(x) (m-x)^2))/length(res))
  }
  rm(res)
  
  # computing observed deviations:
  o <- crossprod(motifMatches,counts)
  expect <- rowMeans(o)
  deviations <- (o-expect)/expect
  
  al <- list( #bg.means=m, bg.sds=sds, o=o, raw.devs=deviations,
    deviations=deviations-m, z=(deviations-m)/sds )
  se <- SummarizedExperiment(al)
  if(!is.null(CD)) colData(se) <- CD
  se
}

.computeDeviationsMemoryUsage <- function(counts, motifMatches, backgrounds,
                                          nthreads=1){
  # simple algorithm memory usage:
  os1 <- 8*ncol(counts)*ncol(motifMatches)*(ncol(backgrounds)+2)/1e+09
  # using Welford's algorithm:
  os2 <- 8*ncol(counts)*ncol(motifMatches)*6*nthreads/1e+09
  list(standard=os1, welford=os2)
}

# computes background deviations using Welford's online SD & mean calculation
.computeBgDevChunkWelford <- function(bgindices, backgrounds, motifMatches,
                                      counts, doVariant=FALSE){
  n <- 0L
  m <- m2 <- sd <- matrix(0, nrow=ncol(motifMatches), ncol=ncol(counts))
  
  for(i in bgindices){
    if(doVariant){
      o <- t(crossprod(counts[backgrounds[,i],],motifMatches))
    }else{
      o <- crossprod(motifMatches[backgrounds[,i],],counts)
    }
    expect <- rowMeans(o)
    o <- (o-expect)/expect
    # Welford's online algorithm:
    n <- n + 1L
    delta1 <- o - m
    m <- m + delta1/n
    delta2 <- o - m
    m2 <- m2 + delta1 * delta2
  }
  list(m, m2/(n-1))
}
