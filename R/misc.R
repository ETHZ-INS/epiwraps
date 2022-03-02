.maketrans <- function(tcol, alpha=100){
  c <- col2rgb(tcol)
  rgb(c["red", 1][[1]], c["green", 1][[1]], c["blue", 1][[1]], 
      alpha, maxColorValue = 255)
}

.parseFiletypeFromName <- function(x, stopOnUnrecognized=TRUE, grOk=FALSE,
                                   trackOk=FALSE){
  if(is(x,"GRanges")){
    if(grOk) return("GRanges")
    stop("The argument should be a filepath!")
  }
  if(inherits(x, "GdObject")){
    if(trackOk) return("track")
    stop("The argument should be a filepath!")
  }
  if(length(x)>1) return(sapply(x, .parseFiletypeFromName))
  if(grepl("\\.bigwig$|\\.bw$", x, ignore.case=TRUE)) return("bw")
  if(grepl("\\.bam$", x, ignore.case=TRUE)) return("bam")
  if(grepl("\\.bed$|\\.narrowpeak$|\\.broadpeak$|\\.gappedpeak$", x, 
           ignore.case=TRUE)) return("bed")
  if(stopOnUnrecognized) stop("Format unrecongized for file\n",x)
  return(NULL)
}

.cleanFileNames <- function(x){
  pat <- paste0("\\.bigwig$|\\.bw$|\\.bam$|\\.bed$|\\.narrowpeak$",
                "|\\.broadpeak$|\\.gappedpeak$")
  gsub(pat, "", basename(x), ignore.case=TRUE)
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


.getBP <- function(x, progressbar=NULL, ...){
  if(is.numeric(x)){
    if(is.null(progressbar)) progressbar <- x > 1
    if(length(x)==1 && isTRUE(x>0)){
      if(x==1) return(BiocParallel::SerialParam(progressbar=progressbar, ...))
      if(.Platform$OS.type=="unix"){
        return(BiocParallel::MulticoreParam(x, progressbar=progressbar, ...))
      }else{
        return(BiocParallel::SnowParam(x, progressbar=progressbar, ...))
      }
    }
  }else if(inherits(x, "BiocParallelParam")){
    return(x)
  }
  stop("BPPARAM should either be a BiocParallelParam object or a positive",
       "integer indicating the number of threads to use.")
}

.comparableMatrices <- function(ml, checkAttributes=FALSE, verbose=TRUE){
  if(!is.list(ml)) ml <- list(signal=ml)
  if(!all(unlist(lapply(ml, FUN=is, class2="normalizedMatrix"))))
    stop("The object is not a list of signal matrices (i.e. normalizedMatrix)")
  dims <- do.call(rbind, lapply(ml, dim))
  if(length(unique(dims[,1]))>1){
    if(length(tt <- table(unlist(lapply(ml, row.names))))==0 ||
       length(i <- names(tt)[which(tt==length(ml))])<2)
      stop("The matrices do not have the same number of regions, and appear to
           have no named region in common.")
    if(verbose)
      warning("The matrices contain different numbers of regions, and only ",
            "common region names (", length(i), ") will be retained.")
    ml <- lapply(ml, FUN=function(x) x[i,])
  }
  if(checkAttributes){
    if(length(unique(dims[,2]))>1)
      stop("The matrices do not have the same numbers of columns/windows!")
    a <- attributes(ml[[1]])
    ai <- c("upsteam_index", "downsteam_index", "target_index", "extend")
    if(!all(unlist(lapply(ml, FUN=function(x){
      identical(attributes(x)[ai],a[ai]) && nrow(x)==nrow(ml[[1]])
    })))) stop("The matrices do not have an homogeneous design!")
  }
  ml
}