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


#' motifFootprint
#' 
#' This is a wrapper around `ATACseqQC::factorFootprints`, making it compatible
#' with a broader variety of inputs and use cases.
#'
#' @param bamfile The path to a bam file
#' @param motif The motif (as probability matrix)
#' @param motif_occurences A GRanges of the motif occurences
#' @param genome Optional genome used for seqlengths (otherwise estimated from
#'   the occurences)
#' @param around How much around the motif to plot
#'
#' @return A plot
#' @export
#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqlengths seqlengths<-
motifFootprint <- function(bamfile, motif, motif_occurences, genome=NULL,
                           around=100){
  stopifnot(is.matrix(motif))
  stopifnot(is(motif_occurences, "GRanges"))
  if(!requireNamespace("ATACseqQC", quietly=TRUE))
    stop("This requires the ATACseqQC package.")
  seqlevelsStyle(motif_occurences) <- seqlevelsStyle(BamFile(bamfile))
  if(all(is.na(seqlengths(motif_occurences)))){
    if(!is.null(genome)){
      seqlengths(motif_occurences) <-
        seqlengths(genome)[seqlevels(motif_occurences)]
    }else{
      seqlengths(motif_occurences) <- sapply(
                                split(end(motif_occurences)+as.integer(around),
                                seqnames(motif_occurences)), max)
    }
  }
  if(is.null(motif_occurences$score)) motif_occurences$score <- 1
  ATACseqQC:::factorFootprints(bamfile, pfm=motif, bindingSites=motif_occurences, 
                              upstream=around, downstream=around,
                              seqlev=seqlevels(motif_occurences))
}
