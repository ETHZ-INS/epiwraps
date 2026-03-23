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
#' @param progress Logical; whether to show a progress bar.
#' @param exclude An optional GRanges of regions for which overlapping reads 
#'   should be excluded.
#' @param ... Passed to `fn`
#'
#' @return A list of whatever `fn` returns
#' @export
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @importFrom rtracklayer path import
#' @importFrom BiocParallel bpnworkers bplapply
#' @importFrom pbapply pblapply
#' @examples
#' # generate dummy regions and save them to a temp file:
#' frags <- tempfile(fileext = ".tsv")
#' d <- data.frame(chr=rep(letters[1:2], each=10), start=rep(100*(1:10),2))
#' d$end <- d$start + 15L
#' write.table(d, frags, col.names=FALSE, row.names=FALSE, sep="\t")
#' # tabix-index it
#' frags <- Rsamtools::bgzip(frags)
#' Rsamtools::indexTabix(frags, format = "bed")
#' # now we can do something chunk-wise, e.g. extract coverage:
#' res <- tabixChrApply(frags, fn=coverage)
#' # aggregate the chunk results into an RleList object:
#' reduceRleLists(res)
tabixChrApply <- function(x, fn, keepSeqLvls=NULL, exclude=NULL, only=NULL,
                          BPPARAM=NULL, progress=TRUE, ...){
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
  
  if(is.null(BPPARAM) || BiocParallel::bpnworkers(BPPARAM)==1){
    if(progress) return(pblapply(slvls, FUN=f2, postfn=fn, ...))
    return(lapply(slvls, FUN=f2, postfn=fn, ...))
  }
  bplapply(slvls, FUN=f2, postfn=fn, ..., BPPARAM=BPPARAM)
}
