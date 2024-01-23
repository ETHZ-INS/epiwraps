#' tabixChrApply
#' 
#' Runs a function on reads/fragments from each chromosomes of a Tabix-indexed 
#' fragment file. This is especially used by other functions to 
#' avoid loading all alignments into memory, or to parallelize reads processing.
#'
#' @param x The path to a tabix-indexed bam file, or a TabixFile object.
#' @param fn The function to be run, the first argument of which should be a
#'   `GRanges`
#' @param keepSeqLvls An optional vector of seqLevels to keep
#' @param BPPARAM A `BiocParallel` parameter object for multithreading. Note 
#'   that if used, memory usage will be high; in this context we recommend a 
#'   high `nChunks`.
#' @param only An optional GRanges of regions for which overlapping reads should
#'   be included. If set, all other reads are discarded.
#' @param only 
#' @param ... Passed to `fn`
#'
#' @return A list of whatever `fn` returns
#' @export
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @importFrom rtracklayer path import
#' @importFrom BiocParallel bpnworkers bplapply
#' @importFrom pbapply pblapply
tabixChrApply <- function(x, fn, keepSeqLvls=NULL, exclude=NULL, only=NULL,
                          BPPARAM=SerialParam(), ...){
  x <- TabixFile(x)
  if(!is.null(exclude)) stopifnot(is(exclude, "GRanges"))
  if(!is.null(only)) stopifnot(is(only, "GRanges"))
  
  f2 <- function(sn, postfn, ...){
    x <- rtracklayer::import(path(x), format="bed",
                             which=GRanges(sn, IRanges(1,5*10^8)))
    if(!is.null(only)){
      .comparableStyles(x, only)
      x <- x[overlapsAny(x, only)]
    }
    if(!is.null(exclude)){
      .comparableStyles(x, exclude)
      x <- x[!overlapsAny(x, exclude)]
    }
    postfn(x, ...)
  }

  slvls <- seqnamesTabix(x)
  if(!is.null(keepSeqLvls)){
    if(length(missingLvls <- setdiff(keepSeqLvls, slvls)>0))
      stop(paste0(
        "Some of the seqLevels specified by `keepSeqLvls` are not in the data.
The first few are:",
head(paste(missingLvls, collapse=", "), 3)))
    slvls <- keepSeqLvls
  }
  
  if(BiocParallel::bpnworkers(BPPARAM)==1){
    return(pblapply(slvls, FUN=f2, postfn=fn, ...))
  }
  bplapply(slvls, FUN=f2, postfn=fn, ..., BPPARAM=BPPARAM)
}
