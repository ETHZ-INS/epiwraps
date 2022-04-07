#' findMotifInstances
#'
#' A wrapper around `TFBSTools` for scanning motif occurence, handling the 
#' coordinate conversion as `memes` does. This is an alternative to the `memes`
#' package, which is more efficient but requires the meme suite to be available.
#' 
#' @param seqs A set of sequences, e.g. `DNAStringSet`, optionally with 
#'   coordinate as names to enable conversion.
#' @param motif A motif, in any format recognized by `universalmotif`
#' @param keepMatchedSeq Logical; whether to keep the matched sequence.
#' @param ... Passed to `TFBSTools::searchSeq`; can for instance be used to set
#'   the number of threads to use, e.g. with `mc.cores=2`
#'   
#' @importFrom Biostrings DNAStringSet
#'
#' @return A `GRanges` object
#' @export
findMotifInstances <- function(seqs, motif, keepMatchedSeq=FALSE, ...){
  seqs <- DNAStringSet(seqs)
  
  if(!is(motif, "PWMatrix")){
    if(!requireNamespace("universalmotif", quietly=TRUE))
      stop("This requires the universalmotif package.")
    motif <- universalmotif:::convert_motifs(motif, "TFBSTools-PWMatrix")
  }
  if(!requireNamespace("TFBSTools", quietly=TRUE))
    stop("This requires the TFBSTools package.")
  x <- suppressWarnings(TFBSTools:::searchSeq(motif, subject=seqs))
  peaks <- strsplit(gsub("-",":",names(seqs)), ":")
  if(all(lengths(peaks)==3)){ # convert relative coordinates to absolute
    chrs <- sapply(peaks,FUN=function(x) x[1])
    offsets <- sapply(peaks,FUN=function(x) as.integer(x[2]))-1L
    i <- rep(seq_along(x),lengths(x))
    return(GRanges(chrs[i], 
                   IRanges(offsets[i]+as.integer(unlist(lapply(x, start))),
                           offsets[i]+as.integer(unlist(lapply(x, end)))),
                   strand=unlist(lapply(x, strand)),
                   score=as.numeric(unlist(lapply(x, FUN=function(x) score(x))))
    ))
  }
  x <- as(x, "GRanges")
  score(x) <- x$absScore
  keepFields <- c("score","relScore")
  if(keepMatchedSeq) keepFields <- c(keepFields, "siteSeqs")
  mcols(x) <- mcols(x)[,intersect(colnames(mcols(x)), keepFields)]
  x
}


