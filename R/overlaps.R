#' regionUpset
#' 
#' A wrapper around \code{\link[UpSetR]{upset}} for GRanges.
#'
#' @param x A named list of genomic ranges (or paths to bed files)
#' @param reference The method for creating the reference windows ('reduce' or
#'   'disjoin'). Alternatively, a `GRanges` object of reference windows.
#' @param returnList Logical; whether to return the list instead of plotting.
#' @param ignore.strand Logical; whether to ignore strands when computing 
#' overlaps (default FALSE). Strand information is ignored if either of the 
#' compared sets of regions is unstranded.
#' @param maxgap Maximum gap between regions to count as an overlap (see 
#'  \code{\link[GenomicRanges]{findOverlaps-methods}}).
#' @param minoverlap Minimum overlap to count as a match (see 
#'  \code{\link[GenomicRanges]{findOverlaps-methods}}).
#' @param ... Further plotting arguments passed to \code{\link[UpSetR]{upset}}.
#'
#' @return
#' @export
#' @importFrom UpSetR upset fromList
#' @importFrom GenomicRanges reduce disjoin 
#' @importFrom IRanges IRanges overlapsAny
#'
#' @examples
#' # random list of GRanges:
#' grl <- lapply(c(A=10,B=20,C=30), FUN=function(x){
#'   GRanges("seq1", IRanges(runif(x,1,1000), width=20))
#' })
#' regionUpset(grl)
regionUpset <- function(x, reference=c("reduce","disjoin"), returnList=FALSE,
                        ignore.strand=FALSE, maxgap=-1L, minoverlap=0L, ...){
  if(is.character(x)) x <- as.list(x)
  if(is.list(x)){
    if(is.null(names(x)) && all(unlist(lapply(x,is.character))) && 
       all(lengths(x)==1))
      names(x) <- unlist(x)
    if(any(lapply(x,FUN=function(x){
      if(is.character(x)){
        if(length(x)!=0) return(TRUE)
        if(!file.exists(x)) return(TRUE)
      }else if(!is(x,"GRanges")){
        return(TRUE)
      }
      FALSE
    }))) stop("`x` should either be i) a named list of GRanges, ii) a 
              GRangesList, or iii) a named list of character vectors, each of
              length 1 and providing the path to a bed file.")
    x <- lapply(x, FUN=function(x){
      if(is.character(x)){
        if(grepl("\\.rds$",x,ignore.case=TRUE)){
          y <- readRDS(x)
        }else{
          y <- importBedlike(x)
        }
        if(!is(y,"GRanges"))
          stop(paste0(y," does not appear to contain genomic ranges."))
        x <- y
      }
      x
    })
    if(is.list(x) && all(unlist(lapply(x,class2="GRanges",is))))
      stopifnot(all(unlist(lapply(x,class2="GRanges",is))))
    x <- as(x, "GRangesList")
  }
  stopifnot(is(x,"GRangesList"))
  if(!is(reference,"GRanges"))
    reference <- switch(match.arg(reference),
                        reduce=reduce(unlist(x),ignore.strand=ignore.strand),
                        disjoin=disjoin(unlist(x),ignore.strand=ignore.strand))
  x <- lapply(x, FUN=function(x)
    which(overlapsAny(reference, x, ignore.strand=ignore.strand, 
                      maxgap=maxgap, minoverlap=minoverlap)))
  if(returnList) return(x)
  UpSetR::upset(UpSetR::fromList(x), ...)
}