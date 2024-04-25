.maketrans <- function(tcol, alpha=100){
  c <- col2rgb(tcol)
  rgb(c["red", 1][[1]], c["green", 1][[1]], c["blue", 1][[1]], 
      alpha, maxColorValue = 255)
}

.parseFiletypeFromName <- function(x, stopOnUnrecognized=TRUE, mustExist=TRUE,
                                   grOk=FALSE, trackOk=FALSE, covOk=FALSE,
                                   requireUnique=FALSE){
  if(is(x,"GRanges")){
    if(grOk) return("GRanges")
    stop("The argument should be a filepath!")
  }
  if(inherits(x, "GdObject")){
    if(trackOk) return("track")
    stop("The argument should be a filepath!")
  }
  if(is(x,"RleList")){
    if(covOk) return("cov")
    stop("The argument should be a filepath!")
  }
  if(length(x)>1){
    y <- vapply(x, stopOnUnrecognized=stopOnUnrecognized, mustExist=FALSE,
                grOk=grOk, trackOk=trackOk, FUN.VALUE=character(1), 
                FUN=.parseFiletypeFromName)
    if(requireUnique && length(y <- unique(y))>1)
      stop("All inputs files should be of the same format!")
    if(mustExist){
      w <- which(y %in% c("bam","bw","bed"))
      if(any(fe <- !file.exists(x[w])))
        stop("The following filepath(s) appear(s) not to exist:\n",
             paste(x[w[fe]], collapse="\n"))
    }
    return(y)
  }
  if(grepl("\\.bigwig$|\\.bw$", x, ignore.case=TRUE)) return("bw")
  if(grepl("\\.bam$", x, ignore.case=TRUE)) return("bam")
  if(grepl("\\.bed$|\\.narrowpeak$|\\.broadpeak$|\\.gappedpeak$", x, 
           ignore.case=TRUE)) return("bed")
  if(mustExist && !file.exists(x))
    stop("The given filepath appears not to exist.")
  if(stopOnUnrecognized) stop("Format unrecongized for file\n",x)
  return(NULL)
}

.cleanFileNames <- function(x){
  if(!is.character(x))
    stop("Cannot create names for non-character input. Please set names.")
  pat <- paste0("\\.bigwig$|\\.bw$|\\.bam$|\\.bed$|\\.narrowpeak$",
                "|\\.broadpeak$|\\.gappedpeak$")
  gsub(pat, "", basename(x), ignore.case=TRUE)
}

.align2cuts <- function(x, size=1L){
  sort(c(resize(x,size,fix="start"), resize(x,size,fix="end")))
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


#' importBedlike
#'
#' Imports a bed-like file as a GRanges object. Uses
#' `rtracklayer` import functions if possible, and falls
#' back onto an import that's not format-committed 
#' otherwise.
#'
#' @param x The path to a bed or bed-like file (can be 
#' gzipped)
#' @param ... passed to \code{\link[data.table]{fread}}
#'
#' @return A `GRanges` object
#' @export
#' @importFrom data.table fread
#' @importFrom rtracklayer import.bed import.bed15
#' @importFrom GenomicRanges GRanges score<-
#' @importFrom S4Vectors mcols mcols<-
#'
#' @examples
#' # example bed file:
#' filepath <- system.file("extdata/example_peaks.bed", 
#'                         package="epiwraps")
#' b <- importBedlike(filepath)
importBedlike <- function(x, ...){
  y <- try(rtracklayer::import.bed(x), silent=TRUE)
  if(is(y,"try-error"))
    y <- try(rtracklayer::import.bed15(x), silent=TRUE)
  if(is(y,"try-error"))
    y <- try(rtracklayer::import.bed(x, format="narrowPeak"), silent=TRUE)
  if(!is(y,"try-error")) return(y)
  y <- as.data.frame(data.table::fread(x, ...))
  strInfo <- NULL
  if(ncol(y)>=6 && all(y[[6]] %in% c(".","*","+","-"))){
    y[[6]] <- gsub(".","*",y[[6]],fixed=TRUE)
    if(!all(y[[6]]==".")) strInfo <- y[[6]]
    y <- y[,-6]
  }
  gr <- GRanges(y[,1], IRanges(y[,2],y[,3]), strand=strInfo)
  if(ncol(y)>=5){
    if(is.numeric(y[[5]])){
      if(!all(y[[5]]==0L)) score(gr) <- y[[5]]
      y <- y[,-5]
    }else if(all(y[,5]==".")){
      y <- y[,-5]
    }
  }
  if(ncol(y)>=4){
    if(length(unique(y[[4]]))==nrow(y)){
      names(gr) <- y[[4]]
      y <- y[,-4]
    }else{
      if(all(y[,4]==".")) y <- y[,-4]
    }
  }
  mcols(gr) <- cbind(mcols(gr),y[,-1:-3,drop=FALSE])
  gr
}


#' Inject (insert) values at positions in a vector
#'
#' @param what The values to inject
#' @param inWhat The vector in which to inject them
#' @param at The positions in the vector *after* which to inject the values.
#'
#' @return A vector of same mode at `inWhat`, and of length equal to the sum of 
#'   the lengths of `what` and `inWhat`.
#' @export
#'
#' @examples
#' inWhat <- 1:10
#' inject(c(21,22,23), inWhat, at=as.integer(c(0,5,10)))
inject <- function(what, inWhat, at){
  stopifnot(is.integer(at))
  if(length(at)==0L) return(inWhat)
  stopifnot(is.vector(inWhat) && is.vector(what))
  if(length(what)==1L) what <- rep(what, length(at))
  stopifnot(length(what)==length(at))
  if(is.factor(inWhat)){
    stopifnot(all(unique(what) %in% levels(inWhat)))
    return(factor(inject(as.integer(factor(what, levels(inWhat))),
                         as.integer(inWhat), at),
                  levels=seq_along(levels(inWhat)), labels=levels(inWhat)))
  }
  stopifnot(min(at)>=0L & max(at)<=(length(inWhat)))
  stopifnot(mode(what)==mode(inWhat))
  # create output vector
  x <- vector(mode=mode(inWhat), length=length(what)+length(inWhat))
  # get index increases
  if(length(at)==1){
    ii <- as.integer(seq_along(inWhat)>at)
  }else{
    ii <- as.integer(cut(seq_along(inWhat), 
                         unique(c(0L,at,length(inWhat))), labels=FALSE))
    ii[is.na(ii)] <- 0L
    if(!any(at==0L)) ii <- ii-1L
  }
  # inject original values
  x[seq_along(inWhat)+ii] <- inWhat
  # inject new values
  x[at+seq_along(at)] <- what
  x
}


.checkMissingSeqLevels <- function(x, seqlvls, fatal=TRUE){
  if(!is.vector(x)) x <- seqlevels(x)
  if(is.integer(x)) x <- names(x)
  if(is.null(seqlvls)) return(invisible(x))
  x <- intersect(x, seqlvls)
  if(length(x)==0) stop("No seqLevel retained!")
  if(length(missingLvls <- setdiff(seqlvls, x))==0) return(invisible(x))
  msg <- paste0(
    "Some of the seqLevels specified by `keepSeqLvls` are not in the data.
The first few are:",
head(paste(missingLvls, collapse=", "), 3))
  if(fatal) stop(msg)
  warning(msg)
  invisible(x)
}


# checks that the regions are compatible with the bw/bam files
.checkRegions <- function(tracks, regions, verbose=TRUE, trimOOR=FALSE){
  if( (is.list(tracks) || is.character(tracks)) && length(tracks)>1 ){
    r2 <- lapply(tracks, regions=regions, verbose=FALSE, trimOOR=trimOOR,
                 FUN=.checkRegions)
    isIn <- do.call(cbind, lapply(r2, FUN=function(x){
      overlapsAny(regions,x,type="equal")
    }))
    r2 <- regions[which(rowSums(isIn)==ncol(isIn))]
  }else{
    if(is.list(tracks)) tracks <- tracks[[1]]
    if(is(tracks, "GRanges")){
      sl <- seqlengths(tracks)
    }else if(.parseFiletypeFromName(tracks)=="bw"){
      sl <- seqlengths(BigWigFile(tracks))
    }else if(.parseFiletypeFromName(tracks)=="bam"){
      sl <- seqlengths(BamFile(tracks))
    }else{
      return(regions)
    }
    if(all(is.na(sl))) return(regions)
    r2 <- keepSeqlevels(regions, intersect(seqlevels(regions), names(sl)),
                        pruning.mode="coarse")
    if(any(is.na(seqlengths(regions))))
      seqlengths(regions) <- sl[seqlevels(regions)]
    if(trimOOR){
      r2b <- trim(r2)
      if(verbose && (length(r2b)!=length(r2) || any(ranges(r2b)!=ranges(r2))))
        message("Some of the regions were out of range and were trimmed.")
      r2 <- r2b
    }
  }
  lost <- length(regions)-length(r2)
  lostp <- round(100*lost/length(regions))
  if(lostp>5 || (verbose && lost>0)){
    msg <- paste0(lost," region(s) (",lostp,"%) were excluded because they ",
                  "were out of range of (some of) the file(s).
This usually happens when the genome annotation used for the files ",
"differs from that on which the regions were based.")
    if(lostp<=5){
      message(msg)
    }else if(lostp<=95){
      warning(msg)
    }else{
      stop(msg)
    }
  }
  return(r2)
}



#' views2Matrix
#'
#' converts a RleViews or RleViewsList with views of the same width to a matrix,
#' setting out-of-bounds regions to NA (or `padVal`).
#'
#' @param v A `RleViews` or `RleViewsList` object with views of the same width.
#' @param padVal The value to assign to out-of-bound regions.
#'
#' @return A numeric matrix.
#' @export
#'
#' @examples
#' # we create an example RleViews with out-of-bound regions:
#' library(IRanges)
#' co <- Rle(values=c(0,1,0), lengths=c(100,50,20))
#' v <- Views(co, c(25,150),c(50,175))
#' # convert to matrix
#' views2Matrix(v)
views2Matrix <- function(v, padVal=NA_integer_){
  if(!is(v, "RleViewsList")) v <- RleViewsList(v)
  ws <- width(v)[[head(which(lengths(v)>0),1)]][[1]]
  stopifnot(all(unlist(width(v))==ws))
  x <- Reduce(c, lapply(v[lengths(v)>0], padVal=padVal, FUN=.view2paddedAL))
  matrix(unlist(x), byrow=TRUE, ncol=ws)
}

# converts a RleViews to an AtomicList, setting out-of-bounds regions to padVal
.view2paddedAL <- function(v, padVal=NA_integer_, forceRetAL=TRUE){
  v2 <- trim(v)
  if(isInt <- is.integer(runValue(v2[[1]]))){
    stopifnot(is.integer(padVal))
  }else{
    stopifnot(is.numeric(padVal))
  }
  toAL <- function(v){
    if(isInt) return(IntegerList(v))
    NumericList(v)
  }
  if(all(width(v2)==width(v))){
    if(forceRetAL) return(toAL(v))
    return(v)
  }
  if(any(w <- width(v2)==0)){
    v <- v[which(!w)]
    v2 <- v2[which(!w)]
    warning(sum(w), " views were excluded as completely out of range.")
  }
  # figure out how much is trimmed on either side
  pleft <- start(v2)-start(v)
  pright <- end(v)-end(v2)
  # concatenate the list elements with their padding
  n <- seq_along(v2)
  v <- splitAsList(c( rep(padVal, sum(pleft)),
                      unlist(toAL(v2), use.names=FALSE),
                      rep(padVal, sum(pright))),
                   c(rep(n, pleft), rep(n, width(v2)), rep(n, pright)))
  names(v) <- names(v2)
  v
}

.comparableStyles <- function(a,b,stopIfNot=TRUE){
  if(.getSeqLevelsStyle(a)==.getSeqLevelsStyle(b)) return(TRUE)
  msg <- paste("It seems your are providing objects for which the seqlevel ",
               "styles do not match.")
  if(stopIfNot) stop(msg)
  warning(msg)
  FALSE
}
.getSeqLevelsStyle <- function(x){
  if(is.character(x)){
    if(grepl("\\.bw$|\\.bigwig", x, ignore.case=TRUE)){
      return(seqlevelsStyle(rtracklayer::BigWigFile(x)))
    }else if(grepl("\\.bam$", x, ignore.case=TRUE)){
      return(seqlevelsStyle(Rsamtools::BamFile(x)))
    }
    return(tryCatch({
      seqlevelsStyle(rtracklayer::BEDFile(x))
    }, error=function(e)
      stop("Unknown filetype.")))
  }
  seqlevelsStyle(x)
}

.filterFrags <- function(frags, only, exclude){
  if(!is.null(only)){
    .comparableStyles(frags, only)
    frags <- frags[overlapsAny(frags, only)]
  }
  if(!is.null(exclude)){
    .comparableStyles(frags, exclude)
    frags <- frags[!overlapsAny(frags, exclude)]
  }
  frags
}

.mdsSortRows <- function(x, scale=FALSE){
  x <- as.matrix(x)
  if(ncol(x)==1) return(order(x))
  if(scale){
    x <- t(x)
    sds <- matrixStats::colSds(x,na.rm=TRUE)
    sds[is.na(sds)] <- 1
    x <- t((x-colMeans(x,na.rm=TRUE))/sds)
  }
  emb <- stats::cmdscale(dist(x), k=2L)
  # taken from package seriation:
  alpha <- atan2(emb[,1], emb[,2])
  o <- order(alpha)
  cut <- which.max(abs(diff(c(alpha[o], alpha[o[1]] + 2 * pi))))
  if(cut != length(o)) o <- o[c((cut + 1):length(o), 1:(cut))]
  o
}