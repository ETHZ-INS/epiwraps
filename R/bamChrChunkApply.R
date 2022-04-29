#' bamChrChunkApply
#' 
#' Runs a function on reads/fragments from chunks of (chromosomes) of an 
#' indexed bam file. This is especially used by other functions to avoid 
#' loading all alignments into memory, or to parallelize reads processing.
#'
#' @param x A bam file.
#' @param FUN The function to be run, the first argument of which should be a
#'   `GRanges`
#' @param paired Logical; whether to consider the reads as paired (fragments, 
#'   rather than reads, will be returned)
#' @param keepSeqLvls An optional vector of seqLevels to keep
#' @param flgs `scanBamFlag` for filtering the reads
#' @param mapqFilter Integer of the minimum mapping quality for reads to be 
#'   included.
#' @param nChunks The number of chunks to use (higher will use less memory but
#'   increase overhead)
#' @param BPPARAM A `BiocParallel` parameter object for multithreading. Note 
#'   that if used, memory usage will be high; in this context we recommend a 
#'   high `nChunks`.
#' @param ... Passed to `FUN`
#'
#' @return A list of whatever `FUN` returns
#' @export
bamChrChunkApply <- function(x, FUN, paired=FALSE, keepSeqLvls=NULL, 
                             flgs=scanBamFlag(), mapqFilter=NA_integer_,
                             nChunks=4, BPPARAM=SerialParam(), ...){
  
  param <- .getBamChunkParams(x, flgs=flgs, keepSeqLvls=keepSeqLvls, 
                              nChunks=nChunks)
  f2 <- function(p, ...){
    if(paired){
      x <- readGAlignmentPairs(x, param=p)
      x <- as(x[isProperPair(x)], "GRanges")
    }else{
      x <- GRanges(readGAlignments(x, param=p))
    }
    FUN(x, ...)
  }
  if(BiocParallel::bpnworkers(BPPARAM)==1){
    return(lapply(param,FUN=f2,...))
  }
  bplapply(param, FUN=f2, ..., BPPARAM=BPPARAM)
}

.getBamChunkParams <- function(x, flgs=scanBamFlag(), keepSeqLvls=NULL, 
                               nChunks=4, ...){
  seqs <- Rsamtools::scanBamHeader(x)[[1]]$targets
  if(!is.null(keepSeqLvls)){
    if(length(missingLvls <- setdiff(keepSeqLvls, names(seqs))>0))
      stop(paste0(
        "Some of the seqLevels specified by `keepSeqLvls` are not in the data.
The first few are:",
head(paste(missingLvls, collapse=", "), 3)))
    seqs <- seqs[keepSeqLvls]
  }
  # generate list of reading params
  if(isTRUE(nChunks) || nChunks>=2){
    if(is.numeric(nChunks) && nChunks>=2){
      stopifnot(round(nChunks)==nChunks)
      nChunks <- as.integer(nChunks)
      seqs <- sort(seqs, dec=TRUE)
      chrGroup <- head(rep(seq_len(nChunks), 
                           ceiling(length(seqs)/nChunks)), length(seqs))
      chrGroup <- split(seqs, chrGroup)
    }else{
      chrGroup <- split(seqs, names(seqs))
    }
    param <- lapply(chrGroup, FUN=function(x)
      ScanBamParam(flag=flgs, ..., which=GRanges(names(x), IRanges(1L,x))))
  }else{
    chrGroup <- list(wg=seqs)
    param <- list(wg=ScanBamParam(flag=flgs, ..., 
                                  which=GRanges(names(seqs), IRanges(1L,seqs))))
  }
  param
}
