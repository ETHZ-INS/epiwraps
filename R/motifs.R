#' motifCoOccurence
#' 
#' Returns regions that have matches for given pairs of motifs within certain 
#' distances of each other.
#'
#' @param motifs A named `PFMatrixList` or `PWMatrixList` object containing 
#'   motifs (only those specified in `pairs` will be used). If you're not 
#'   familiar with these objects, see the `TFBSTools` package, and the 
#'   `univervalmotif` package on how to convert motifs.
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
#'   quantiles for each pair of motifs. Note that `minDist` is ignored in this
#'   mode.
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
#' @importFrom Rsamtools FaFile
motifCoOccurence <- function(motifs, pairs, regions, genome, centerDist=TRUE,
                             minDist=5, maxDist=50, exclusiveDist=TRUE,
                             restrictToRegions=FALSE,
                             nDistQuantiles=NULL, ignore.strand=TRUE){
  
  stopifnot(is(motifs, "XMatrixList"))
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
  if(!all(sapply(pairs, \(x) x %in% names(motifs)))){
    stop("Some of the motifs specific in `pairs` appear not to be in `motifs`.")
  }
  if(is.null(names(pairs)))
    pairs <- setNames(pairs, sapply(pairs, paste, collapse="+"))
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


#' getInteractionsDeviations
#'
#' Get \code{\link[betterChromVAR]{betterChromVAR}} deviations for overlaps 
#' between motifs and a bait motif.
#' 
#' @param se An object inheriting `SummarizedExperiment`, containing a `counts`
#'   assay.
#' @param annotation A (sparse) matrix or object inheriting 
#'   `SummarizedExperiment` containing motif annotations for each row of `se`.
#' @param bait The name of the bait motif (corresponding to a column of 
#'   `annotation`), for which to search for interacting motifs.
#' @param minCount The minimum overlap size for a motif pair to be considered.
#' @param ... Passed to `betterChromVAR`
#'
#' @returns A `SummarizedExperiment` of deviations for the intersection of 
#'   all motifs with the bait motif.
#' @importFrom betterChromVAR betterChromVAR
#' @export
getInteractionsDeviations <- function(se, annotation, bait, minCount=20, ...){
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(nrow(se)==nrow(annotation))
  stopifnot(length(bait)==1 && bait %in% colnames(annotation))
  if(inherits(annotation, "SummarizedExperiment")){
    annotation <- assay(annotation)
  }
  moi2 <- annotation & annotation[,bait]
  moi2 <- moi2[,colSums(moi2)>=minCount]
  ov <- colSums(moi2)
  dev <- betterChromVAR(se, moi2, ...)
  tot <- colSums(annotation[,colnames(moi2)])
  totBait <- as.integer(tot[bait])
  expected <- totBait*tot/nrow(annotation)
  rowData(dev)$overlap.percentOfBait <- round(100*(ov/totBait),2)
  jac <- colOverlaps(annotation[,colnames(moi2)], annotation[,bait],
                     metric="jaccard")
  rowData(dev)$overlap.jaccard <- as.numeric(jac)
  rowData(dev)$overlap.log2Enr <- log2(ov/expected)
  rowData(dev)$isBait <- row.names(dev)==bait
  dev
}

#' discoverMotifInteractions
#' 
#' Discover motifs interacting with a bait motif, based on the betterChromVAR
#' package.
#'
#' @param dev A deviations SummarizedExperiment produced by 
#'   \code{\link{getInteractionsDeviations}}
#' @param group The column of `colData(dev)` containing the experimental groups
#'   of interest. To column will be coerced to a factor, and all comparisons 
#'   will be made relative to the first level of the factor.
#' @param covar Optional columns of `colData(dev)` containing covariates to 
#'   correct for.
#' @param global Logical; whether to test the contrasts globally, rather than
#'   individually. Has no effect when there is a single contrast 
#'   (i.e. two-group comparison).
#' @param weights The weights to use, either 'fixed' (uses the sqrt number of 
#'   sites), 'trend' (default; uses the trend over the number of sites), or 
#'   'none' (no weights; not recommended). We recommend 'trend' if the number
#'   of motifs is sufficiently large (e.g. >100).
#' @param useAssay Which assay to use. We recommend setting 'adjZ' (default),
#'   which uses z-scores but corrects the scale of the bait's z-scores to its
#'   expectation for the intersection size.
#' @author Pierre-Luc Germain
#'
#' @returns A data.frame.
#' @importFrom limma lmFit makeContrasts eBayes topTable contrasts.fit
#' @importFrom stats model.matrix
#' @export
discoverMotifInteractions <- function(dev, group, covar=c(), global=FALSE,
                                      weights=c("fixed","trend","none"),
                                      useAssay=c("adjZ","deviations","z")){
  stopifnot(inherits(dev, "SummarizedExperiment"))
  stopifnot(sum(rowData(dev)$isBait)==1)
  useAssay <- match.arg(useAssay)
  weights <- match.arg(weights)
  bait <- row.names(dev)[which(rowData(dev)$isBait)]
  wBait <- which(row.names(dev)==bait)
  N <- rowData(dev)$N
  if(useAssay=="adjZ"){
    useAssay <- "z"
    # shrink the bait's z to expectation for the interaction size
    nf <- sqrt(N/N[wBait])
    e1 <- outer(nf[-wBait], assay(dev,useAssay)[bait,])
  }else{
    e1 <- matrix(assay(dev,useAssay)[bait,], nrow=nrow(dev)-1, ncol=ncol(dev))
  }
  e2 <- assay(dev,useAssay)[-wBait,]
  cd <- as.data.frame(rbind(colData(dev),colData(dev)))
  cd$isInt <- rep(c(FALSE,TRUE), each=ncol(dev))
  cd$combined <- factor(paste(cd[[group]], cd$isInt, sep="_"))
  cd[[group]] <- factor(cd[[group]])
  formula <- paste(c("~0+combined", covar), collapse="+")
  mm <- model.matrix(as.formula(formula), data=cd)
  colnames(mm) <- gsub("combined","",colnames(mm))
  e <- cbind(e1,e2)
  w <- NULL
  if(weights=="fixed"){
    motif_weight <- sqrt(rowData(dev)$N)
    w <- matrix(motif_weight[-wBait], nrow=nrow(e), ncol=ncol(e))
    w[, seq_len(ncol(dev))] <- motif_weight[wBait]
    fit <- lmFit(e, mm, weights=w)
  }else{
    fit <- lmFit(e, mm)
  }
  refG <- levels(cd[[group]])[1]
  grs <- levels(cd[[group]])[-1]
  names(grs) <- paste0(grs,"-",refG)
  contrast_strings <- sapply(grs, function(g) {
    paste0("(", g, "_TRUE - ", g, "_FALSE) - (", 
           refG, "_TRUE - ", refG, "_FALSE)")
  })
  cont_matrix <- makeContrasts(contrasts=contrast_strings, levels=mm)
  if(weights=="trend"){
    motif_weight <- sqrt(rowData(dev)$N[-wBait])
    fit <- eBayes(contrasts.fit(fit, cont_matrix), trend=motif_weight)
  }else{
    fit <- eBayes(contrasts.fit(fit, cont_matrix))
  }
  
  if(!isTRUE(global)){
    res <- dplyr::bind_rows(lapply(setNames(names(grs),names(grs)), \(co){
      res <- as.data.frame(topTable(fit, coef=co, number=Inf))
      colnames(res)[1] <- "difference"
      cbind(motif=row.names(res), res)
    }), .id="contrast")
  }else{
    res <- as.data.frame(topTable(fit, coef=names(grs), number=Inf))
    res <- cbind(motif=row.names(res), res)
  }
  res$B <- res$AveExpr <- NULL
  cols <- c("N",grep("^overlap", colnames(rowData(dev)), value=TRUE))
  d <- as.data.frame(rowData(dev)[,cols])
  colnames(d)[1] <- "overlap.N"
  cbind(res[,c("contrast","motif")], d[res$motif,],
        res[,setdiff(colnames(res), c("contrast", "motif"))])
}

#' exploreMotifInteraction
#' 
#' We compare deviations in sites that harbor a single of the two motifs, or
#' the pair of motifs at given distances from each other.
#'
#' @param counts An object inheriting from `RangedSummarizedExperiment`, with
#'   a `bias` column in its `rowData`.
#' @param motifs A named pair of motifs, in a format recognized by 
#'   \code{\link[motifmatchr]{matchMotif}}.
#' @param genome A BSgenome or FaFile object.
#' @param nDistQuantiles The number of distance quantiles between the motif 
#'   matches (up to `maxDist`)
#' @param maxDist The maximum distance between matches.
#' @param ... Any further argument passed to 
#'   \code{\link[betterChromVAR]{betterChromVAR}}.
#' @author Pierre-Luc Germain
#'
#' @returns A SummarizedExperiment object of deviations for different types of
#'   interactions between the motifs.
#' 
#' @details
#' Motif pairs are split according to the distance between their centers up to
#' `maxDist` (default 300) into `nDistQuantiles` (default 4) similarly-populated
#' quantiles. Deviations are reported for these different bins, as well as sites
#' harboring only a single motif.
#' @export
#' @importFrom betterChromVAR computeDeviationsAnalytic normalizeDevsForSize
#' @importFrom motifmatchr matchMotifs
#' @examples
exploreMotifInteraction <- function(counts, motifs, genome, bg=NULL,
                                    nDistQuantiles=5L, maxDist=500L, ...){
  stopifnot(length(motifs)==2)
  stopifnot(inherits(counts, "RangedSummarizedExperiment"))
  moi <- motifmatchr::matchMotifs(motifs, counts, genome=genome)
  moi2 <- motifCoOccurence(motifs, list(names(motifs)), rowRanges(counts),
                           genome = genome, nDistQuantiles=nDistQuantiles,
                           maxDist=maxDist)
  breaks <- round(attr(moi2, "breaks")[[1]][-1])
  bins <- c(paste0("d<=", breaks[1]), paste0(breaks[-length(breaks)],
                                             "<d<=", breaks[-1]))
  moi2flat <- Reduce(cbind, moi2)
  colnames(moi2flat) <- paste(colnames(moi2[[1]]),
                              rep(names(moi2), each=ncol(moi2[[1]])), sep=".")
  moiSinglet <- (assay(moi)-(Reduce("+", lapply(moi2,Matrix::rowSums))>0))>0
  colnames(moiSinglet) <- paste0(colnames(moiSinglet),"+,",
                                 rev(colnames(moiSinglet)),"-")
  moi <- cbind(assay(moi), moiSinglet, moi2flat)
  for(b in seq_along(bins)){
    colnames(moi) <- gsub(paste0("\\.bin",b), paste0("\n", bins[b]),
                          colnames(moi))
  }
  if(is.null(bg)){
    dev <- betterChromVAR(counts, moi, ...)
  }else{
    dev <- computeDeviationsAnalytic(counts, bg, moi, ...)
  }
  normalizeDevsForSize(dev)
}
