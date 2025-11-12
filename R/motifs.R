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


#' motifCoOccurence
#' 
#' Returns regions that have matches for given pairs of motifs within certain 
#' distances of each other.
#'
#' @param motifs A list of motifs (only those specified in `pairs` will be 
#'   used).
#' @param pairs A list of pairs of motifs for which to compute co-occurence.
#'   Specifically, this should be a (optionally named) list of character vectors
#'   of length 2, corresponding to names in `motifs`.
#' @param regions The regions in which to search for motif matches.
#' @param genome A genome object or path to a genome fasta file.
#' @param centerDist Logical; whether to compute distances from the center of 
#'   the motifs.
#' @param minDist The minimum distance between matches for a co-occurrence. This
#'   is primarily used to exclude interactions between highly-similar, 
#'   overlapping motifs, and to focus on putative cooperative interactions 
#'   between TFs. If interested in competing TFs, set this to <=0.
#' @param maxDist The maximum distance between the matches for a co-occurence.
#' @param exclusiveDist Logical; whether the co-occurences should be counted 
#'  only for the smallest of the distance thresholds. Ignored if `minDist` and
#'  `maxDist` have a length of 1.
#' @param ignore.strand Logical; whether to ignore the strand for co-occurence
#'   (default TRUE).
#' 
#' @details
#' Note that both `minDist` and `maxDist`, rather than specifying a single 
#' threshold, can compute co-occurence for multiple thresholds (which is must
#' faster than running the function multiple times). `minDist` and `maxDist` 
#' should have the same length, and the corresponding entries will be used 
#' together to produce multiple matrices. (If `exclusiveDist=TRUE`, the 
#' distance threshold are exclusive, meaning that the higher threshold does not
#' include co-occurences of the lower thresholds.)
#' Also note that all matches are stored in memory, so using this function 
#' across the entire genome is not advisable (unless for very few motifs).
#' 
#'
#' @returns A list of sparse logical matrices, with one matrix for each value 
#'   of `minDist`/`maxDist`.
#' @export
#' @importFrom motifmatchr matchMotifs
#' @importFrom TFBSTools PFMatrixList
#' @importFrom universalmotif convert_motifs
#' @importFrom Rsamtools FaFile
motifCoOccurence <- function(motifs, pairs, regions, genome, centerDist=TRUE,
                             minDist=5, maxDist=50, exclusiveDist=TRUE,
                             ignore.strand=TRUE){
  stopifnot(is.list(pairs) && all(lengths(pairs)==2))
  stopifnot(is(regions, "GRanges"))
  stopifnot(length(minDist)==length(maxDist))
  if(!is(motifs, "XMatrixList")){
    motifs <- do.call(TFBSTools::PFMatrixList,
                      convert_motifs(motifs, class="TFBSTools-PFMatrix"))
  }
  if(!all(sapply(pairs, \(x) x %in% names(motifs)))){
    stop("Some of the motifs specific in `pairs` appear not to be in `motifs`.")
  }
  if(is.null(names(pairs))) pairs <- setNames(pairs, sapply(pairs, paste, collapse="+"))
  if(is.character(genome) && length(genome)==1)
    genome <- Rsamtools::FaFile(genome)
  
  matches <- matchMotifs(motifs[unique(unlist(pairs))], regions,
                         genome=genome, out="positions")
  if(centerDist) matches <- lapply(matches, resize, fix="center", width=1)
  distCrit <- setNames(seq_along(minDist), paste0(minDist,"<= d <=",maxDist))
  lapply(distCrit, \(i){
    as(sapply(pairs, \(x){
      o <- findOverlapPairs(matches[[x[1]]], matches[[x[2]]],
                            maxgap=maxDist[i]-1L, ignore.strand=ignore.strand)
      if(exclusiveDist)
        o <- o[which(abs(start(first(o))-start(second(o)))>=minDist[i])]
      gr <- GRanges(seqnames(first(o)),
                    IRanges(start=pmin(start(first(o)), start(second(o))),
                            end=pmax(start(first(o)), start(second(o)))))
      return(overlapsAny(regions, gr))
    }), "sparseMatrix")
  })
}
