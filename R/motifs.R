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
#'   between TFs. If interested in competing TFs, set this to <=0. Note that
#'   `minDist` is ignored when using `nDistQuantiles`.
#' @param maxDist The maximum distance between the matches for a co-occurence.
#' @param exclusiveDist Logical; whether the co-occurences should be counted 
#'  only for the smallest of the distance thresholds. Ignored if `minDist` and
#'  `maxDist` have a length of 1.
#' @param restrictToRegions Logical; whether to restrict the search for motifs 
#'  to `regions`. If FALSE (default), one of the motifs still needs to be within
#'  the regions, but the other can be +/- `maxDist` around the region.
#' @param nDistQuantiles The number of distance quantiles to use. Disabled by
#'   default and requires `exclusiveDist=TRUE`. When a positive integer, the 
#'   span from 0 to `maxDist` is split into `nDistQuantiles` number of 
#'   quantiles for each pair of motifs.
#' @param ignore.strand Logical; whether to ignore the strand for co-occurence
#'   (default TRUE).
#' 
#' @details
#' Note that both `minDist` and `maxDist`, rather than specifying a single 
#' threshold, can compute co-occurence for multiple thresholds (which is must
#' faster than running the function multiple times). `minDist` and `maxDist` 
#' should have the same length, and the corresponding entries will be used 
#' together to produce multiple matrices.
#' If `exclusiveDist=TRUE`, the distance threshold are exclusive, meaning that 
#' the higher threshold does not include co-occurences of the lower thresholds.
#' When using multiple exclusive distance bins, it is also possible to have the
#' bin boundaries adjusted for each pair of motifs so that each been is equally
#' populated by setting `nDistQuantiles` (the breaks used can be retrieved from
#' `attr(results, "breaks")` ).
#' Also note that all matches are stored in memory, so using this function 
#' across the entire genome is not advisable (unless for very few motifs).
#' 
#'
#' @returns A list of sparse logical matrices, with one matrix for each value 
#'   of `minDist`/`maxDist` (or each quantile bin).
#' @export
#' @importFrom motifmatchr matchMotifs
#' @importFrom TFBSTools PFMatrixList
#' @importFrom universalmotif convert_motifs
#' @importFrom Rsamtools FaFile
motifCoOccurence <- function(motifs, pairs, regions, genome, centerDist=TRUE,
                             minDist=5, maxDist=50, exclusiveDist=TRUE,
                             restrictToRegions=FALSE,
                             nDistQuantiles=NULL, ignore.strand=TRUE){
  
  stopifnot(is.list(pairs) && all(lengths(pairs)==2))
  stopifnot(is(regions, "GRanges"))
  stopifnot(length(minDist)==length(maxDist))
  stopifnot(is.numeric(maxDist) && all(maxDist>1))
  stopifnot(is.numeric(minDist) && all(minDist>0))
  
  if(!is.null(nDistQuantiles)){
    stopifnot(is.numeric(nDistQuantiles) && length(nDistQuantiles)==1 &&
                nDistQuantiles>1 && (nDistQuantiles %%1==0))
    if(!exclusiveDist)
      stop("`nDistQuantiles` only makes sense with `exclusiveDist=TRUE`")
    if(length(maxDist)>1)
      stop("`maxDist` and `minDist` should each have a length of 1 when using",
           " `nDistQuantiles`.")
  }
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
  
  if(restrictToRegions){
    r2 <- regions
  }else{
    r2 <- resize(regions, width=width(regions)+2*max(maxDist), fix="center")
  }
  matches <- matchMotifs(motifs[unique(unlist(pairs))], r2,
                         genome=genome, out="positions")
  
  if(centerDist) matches <- lapply(matches, resize, fix="center", width=1)
  
  if(is.null(nDistQuantiles)){
    distCrit <- setNames(seq_along(minDist), paste0(minDist,"<= d <=",maxDist))
    return(lapply(distCrit, \(i){
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
    }))
  }

  # Using adaptative distance thresholds
  res <- lapply(pairs, \(x){
    
    out <- list(bins=vector(mode="integer", length=length(regions)), q=NULL)
    
    # for each region, find the shortest distance between the motifs
    # we need to compute shortest distances in both directions:
    fn <- function(m1,m2){
      d1 <- as.data.frame(distanceToNearest(m1,m2,ignore.strand=ignore.strand))
      d1 <- d1[which(abs(d1$distance)<=(maxDist)),]
      # map back to regions:
      o <- findOverlaps(m1, regions, ignore.strand=ignore.strand)
      o <- o[!duplicated(from(o))]
      m1$region <- NA_integer_
      m1$region[from(o)] <- to(o)
      d1$peak <- m1$region[d1[,1]]
      d1
    }
    d1 <- fn(matches[[x[1]]], matches[[x[2]]])
    d2 <- fn(matches[[x[2]]], matches[[x[1]]])
    # merge the two, and find the shortest distance for each region
    colnames(d1) <- c("r1","r2","d1","peak1")
    colnames(d2) <- c("r2","r1","d2","peak2")
    m <- merge(d1, d2, by=c("r1","r2"), all=TRUE)
    if(nrow(m)<2) return(out)
    m$peak <- ifelse( is.na(m$peak1) | (!is.na(m$peak2) & 
                                          !isFALSE(abs(m$d1)>abs(m$d2))),
                      m$peak2, m$peak1)
    m$distance <- pmin(m$d1,m$d2,na.rm=TRUE)
    m <- m[which(!is.na(m$peak) & !is.na(m$distance)),]
    ag <- aggregate(m$distance, by=list(peak=m$peak), FUN=min, na.rm=TRUE)
    if(nrow(ag)<2) return(out)
    
    # define the distance quantiles
    q <- unique(quantile(ag$x, c(0,seq_len(nDistQuantiles))/nDistQuantiles))
    if(length(q)==1) return(out)
    ag$bin <- cut(ag$x, breaks=q, labels=FALSE, include.lowest=TRUE)
    out$bins[ag$peak] <- ag$bin
    out$q <- as.numeric(q)
    out
  })
  res <- res[!sapply(res, \(x) is.null(x$q))]
  bins <- seq_len(nDistQuantiles);
  names(bins) <- paste0("bin", bins)
  a <- lapply(bins, \(i){
    as(sapply(res, \(x) x$bins==i), "sparseMatrix")
  })
  attr(a, "breaks") <- lapply(res, \(x) x$q)
  a
}
