#' getCovStats
#'
#' @param x A (named) vector of paths to bigwig files
#' @param binSize The size of bins
#' @param nbBins The number of random bins. More bins gives more accurate 
#' readouts but take longer to compute.
#' @param exclude Region to exclude
#' @param genome A named vector of chromosome sizes. If a GRanges without 
#' seqlengths, the maximum coordinate of each chromosome will be used.
#' @param BPPARAM BioParallel BPPARAM for multithreading across files.
#'
#' @return A named list of QC tables
#' @export
#' @importFrom BiocParallel bplapply SerialParam
#' @import GenomicRanges S4Vectors
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom dplyr bind_rows
#' @importFrom rtracklayer import BigWigSelection
getFingerprints <- function(x, binSize=1000, nbBins=10000, exclude=NULL, 
                            genome=NULL, BPPARAM=SerialParam()){
  if(is(genome, "GRanges")){
    maxes <- seqlengths(genome)
    if(any(is.na(maxes)))
      maxes <- vapply(split(end(genome), 
                            droplevels(as.factor(seqnames(genome)))),
                      FUN=max, FUN.VALUE=integer(1L))
  }else if(is(genome, "EnsDb")){
    maxes <- genome(genome)
  }else{
    maxes <- genome
  }
  if(is.null(names(x)))
    names(x) <- make.unique(gsub("\\.bw|\\.bigwig", "", basename(x), 
                                 ignore.case=TRUE))
  chr <- table(sample(factor(names(maxes),names(maxes)), size=nbBins, 
                      prob=maxes/sum(maxes), replace=TRUE))
  pos <- lapply(seq_along(chr), FUN=function(x)
                      sample.int(maxes[x]-binSize,chr[[x]], replace=TRUE))
  gr <- sort(GRanges(rep(names(chr),lengths(pos)), 
                     IRanges(unlist(pos), width=binSize)))
  covs <- bplapply(x, BPPARAM=BPPARAM, FUN=function(x){
    x <- rtracklayer::import(x, format="BigWig", selection=BigWigSelection(gr))
    o <- findOverlaps(gr,x)
    tmp <- max(splitAsList(x$score[to(o)], from(o)))
    if(length(tmp)!=length(gr)){
      y <- rep(0,length(gr))
      y[as.integer(names(tmp))] <- as.numeric(tmp)
    }else{
      y <- as.numeric(tmp)
    }
    y
  })
  if(length(covs)>1){
    tmp <- do.call(cbind, covs)
    pea <- cor(tmp)
    spea <- cor(tmp, method="spearman")
    rm(tmp)
  }
  covs <- lapply(covs, sort)
  d1 <- dplyr::bind_rows(lapply(covs, FUN=function(x){
    data.frame(coverage=unique(x),
               fraction.above=vapply(unique(x), FUN.VALUE=numeric(1),
                                     FUN=function(y) sum(x>y)/length(x)))
  }), .id="file")
  d2 <- dplyr::bind_rows(lapply(covs, FUN=function(x){
    data.frame(rank=rank(x), fraction.of.highest=x/max(x))
  }), .id="file")
  ll <- list(coverage=d1, enrichment=d2)
  if(length(covs)>1) ll <- c(ll, list(cor.pearson=pea, cor.spearman=spea))
  ll
}

#' plotCovStats
#' 
#' Plots coverage statistics, such as as fingerprint plot.
#'
#' @param qc A list of coverage statistics, as produced by 
#'  \code{\link{getCovStats}}.
#' @param labels Passed to \code{\link[cowplot]{plot_grid}}.
#'
#' @return A grid object to be plotted.
#' @export
#' @import ggplot2
#' @importFrom cowplot plot_grid get_legend
plotCovStats <- function(qc, labels="AUTO", show.legend=TRUE){
  p1 <- ggplot(qc$coverage, aes(coverage, fraction.above, colour=file)) + 
    geom_line() + labs(x="Read density", y="Fraction of regions > density") +
    theme(legend.position="top", legend.direction = "horizontal")
  qc$enrichment$rank <- qc$enrichment$rank/max(qc$enrichment$rank)
  p2 <- ggplot(qc$enrichment, aes(rank, fraction.of.highest, colour=file)) + 
    geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey") +
    geom_line() + labs(x="Region rank", y="Fraction of highest density") +
    theme(legend.position="none")
  p <- plot_grid(p1 + theme(legend.position="none"), p2, 
                 labels=labels, nrow=1, scale=0.95)
  if(!show.legend || length(unique(qc$coverage$file))==1) return(p)
  plot_grid(p, get_legend(p1), nrow=2, rel_heights = c(3,1))
}

#' plotCorFromCovStats
#'
#' @param qc A list of coverage statistics, as produced by 
#'   \code{\link{getCovStats}}.
#' @param method The correlation metrics to include
#' @param col Optional heatmap colors
#' @param na_col Color for the diagonal; passed to 
#'   \code{\link[ComplexHeatmap]{Heatmap}}.
#' @param column_title Column title (if NULL, uses the metric)
#' @param ... Passed to \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return A `Heatmap` or `HeatmapList` object ready to be plotted.
#' @export
#' @import ComplexHeatmap
#' @importFrom viridisLite viridis inferno
plotCorFromCovStats <- function(qc, method=c("pearson","spearman"), col=NULL, 
                                na_col="white", column_title=NULL, ...){
  if(is.null(qc$cor.pearson) || is.null(qc$cor.spearman))
    stop("The object does not contain the necessary information - was 
         getCovStats run on multiple samples?")
  method <- match.arg(method, several.ok=TRUE)
  dat <- lapply(setNames(method, method), FUN=function(x){
    m <- switch(x, 
                pearson=qc$cor.pearson,
                spearman=qc$cor.spearman)
    diag(m) <- NA
    m
  })
  o <- as.dendrogram(hclust(dist(do.call(cbind, dat))))
  h <- lapply(names(dat), FUN=function(x){
    if(is.null(col)){
      col = switch(x,
                   pearson=viridisLite::viridis(50),
                   spearman=viridisLite::inferno(50))
    }
    if(is.null(column_title)) column_title <- x
    Heatmap(dat[[x]], name=x, col=col, column_title=column_title,
            cluster_rows=o, cluster_columns=o, na_col="white", ...)
  })
  if(length(h)==1) return(h[[1]])
  hl <- NULL
  for(hi in h) hl <- hl + hi
  hl
}