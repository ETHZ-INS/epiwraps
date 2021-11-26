.maketrans <- function(tcol, alpha=100){
  c <- col2rgb(tcol)
  rgb(c["red", 1][[1]], c["green", 1][[1]], c["blue", 1][[1]], 
      alpha, maxColorValue = 255)
}

.parseFiletypeFromName <- function(x, stopOnUnrecognized=TRUE){
  if(length(x)>1) return(sapply(x, .parseFiletypeFromName))
  if(grepl("\\.bigwig$|\\.bw$|\\.bam$", x, ignore.case=TRUE)) return("bw")
  if(grepl("\\.bam$", x, ignore.case=TRUE)) return("bam")
  if(grepl("\\.bed$|\\.narropeak$|\\.broadpeak$|\\.gappedpeak", x, 
           ignore.case=TRUE)) return("bed")
  if(stopOnUnrecognized) stop("Format unrecongized for file\n",x)
  return(NULL)
}

.align2cuts <- function(x){
  bstart <- bend <- GRanges(x)
  end(bstart) <- start(bstart)
  start(bend) <- end(bend)
  sort(c(bstart,bend))
}


#' breakStrings
#'
#' breaks a string of long-enough words (or vector thereof) into two lines.
#'
#' @param x A character (or factor) vector.
#' @param minSizeForBreak The minimum number of characters to be broken into 
#' two lines.
#' @param lb The separation character (e.g. line break).
#'
#' @return A character vector of length=length(x).
#' @export
#' @examples
#' breakStrings("this is too long for practical purposes")
breakStrings <- function (x, minSizeForBreak=20, lb="\n"){
  if(is.factor(x)){
    levels(x) <- breakStrings(levels(x), minSizeForBreak, lb)
    return(x)
  }
  stopifnot(is.character(x))
  sapply(x, minSizeForBreak = minSizeForBreak, lb = lb,
         FUN=function(x, minSizeForBreak, lb){
    if (nchar(x) <= minSizeForBreak) 
      return(x)
    g <- gregexpr(" ", x)[[1]]
    if (length(g) == 0) 
      return(x)
    if (length(g) == 1 & all(g == -1)) 
      return(x)
    mid <- nchar(x)/2
    mid <- g[order(abs(g - mid))[1]]
    substr(x, mid, mid) <- lb
    return(x)
  })
}