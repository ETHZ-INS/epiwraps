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
