#' tabixChrApply
#' 
#' Runs a function on reads/fragments from each chromosomes of a Tabix-indexed 
#' fragment file. This is especially used by other functions to 
#' avoid loading all alignments into memory, or to parallelize reads processing.
#'
#' @param x The path to a tabix-indexed bam file, or a TabixFile object.
#' @param FUN The function to be run, the first argument of which should be a
#'   `GRanges`
#' @param keepSeqLvls An optional vector of seqLevels to keep
#' @param BPPARAM A `BiocParallel` parameter object for multithreading. Note 
#'   that if used, memory usage will be high; in this context we recommend a 
#'   high `nChunks`.
#' @param exclude An optional GRanges of regions for which overlapping reads 
#'   should be excluded.
#' @param ... Passed to `FUN`
#'
#' @return A list of whatever `FUN` returns
#' @export
#' @importFrom Rsamtools TabixFile seqnamesTabix import
tabixChrApply <- function(x, FUN, keepSeqLvls=NULL, exclude=NULL,
                          BPPARAM=SerialParam(), ...){
  if(!is(x, "TabixFile")) x <- TabixFile(x)
  if(!is.null(exclude)) stopifnot(is(exclude, "GRanges"))
  
  f2 <- function(sn, ...){
    x <- rtracklayer::import(x, format="bed", 
                             which=GRanges(sn, IRanges(1,5*10^8)))
    if(!is.null(exclude)) x <- x[!overlapsAny(x, exclude)]
    FUN(x, ...)
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
    return(lapply(slvls, FUN=f2, ...))
  }
  bplapply(slvls, FUN=f2, ..., BPPARAM=BPPARAM)
}
