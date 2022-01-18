#' getCovStats
#' 
#' Assembles read distribution statistics from a set of bigwig files based on
#' random windows.
#'
#' @param x A (named) vector of paths to bigwig files (all from the same genome)
#' @param binSize The size of bins
#' @param nbBins The approximate number of random bins. More bins gives more 
#' accurate readouts but take longer to read and compute.
#' @param exclude Region to exclude
#' @param canonical.chr Logical; whether to restrict the sampling to standard
#' chromosomes.
#' @param maxCovQuant The quantile to use as maximum coverage (default 0.999)
#' @param BPPARAM BioParallel BPPARAM for multithreading across files.
#'
#' @return A named list of QC tables
#' @export
#' @importFrom BiocParallel bplapply SerialParam
#' @import GenomicRanges S4Vectors
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom dplyr bind_rows
#' @importFrom rtracklayer import BigWigSelection BigWigFile
#' @importFrom Rsamtools BamFile scanBamFlag ScanBamParam
#' 
#' @examples 
#' # we use an example bigwig file
#' bwf <- system.file("extdata/example_atac.bw", package="epiwraps")
#' # because most of the file is empty, we'll exclude some of the ranges
#' cs <- getCovStats(bwf, exclude=GRanges("1", IRanges(1, 4300000)))
#' plotCovStats(cs)
getCovStats <- function(x, binSize=1000, nbBins=10000, exclude=NULL, 
                        canonical.chr=TRUE, maxCovQuant=0.999, 
                        BPPARAM=SerialParam()){
  if(.parseFiletypeFromName(x[[1]])=="bam"){
    maxes <- seqlengths(BamFile(x[[1]]))
  }else{
    maxes <- seqlengths(BigWigFile(x[[1]]))
  }
  stopifnot(length(maxes)>0)
  if(canonical.chr) maxes <- maxes[grep("^Y$|^X$|^[0-9]$", ignore.case=TRUE,
                                        gsub("chr|CHR","",names(maxes)))]
  if(length(maxes)==0)
    stop("No chromosome left after filtering! Consider using ",
         "canonical.chr=FALSE")
  if(is.null(names(x)))
    names(x) <- make.unique(gsub("\\.bw$|\\.bigwig$|\\.bam$", "", basename(x), 
                                 ignore.case=TRUE))
  chr <- table(sample(factor(names(maxes),names(maxes)), size=nbBins, 
                      prob=maxes/sum(maxes), replace=TRUE))
  pos <- lapply(seq_along(chr), FUN=function(x)
                      sample.int(maxes[x]-binSize,chr[[x]], replace=TRUE))
  gr <- sort(GRanges(rep(names(chr),lengths(pos)), 
                     IRanges(unlist(pos), width=binSize)))
  if(!is.null(exclude)) gr <- gr[!overlapsAny(gr, exclude)]
  covs <- bplapply(x, BPPARAM=BPPARAM, FUN=function(x){
    ftype <- .parseFiletypeFromName(x)
    if(ftype=="bam"){
      flgs <- scanBamFlag(isDuplicate=FALSE, isSecondaryAlignment=FALSE)
      params <- Rsamtools::ScanBamParam(which=gr, flag=flgs)
      x <- as(GenomicAlignments::readGAlignmentPairs(x, param=params), "GRanges")
      y <- as.integer(countOverlaps(gr,x))
    }else if(ftype=="bw"){
      x <- rtracklayer::import(x, format="BigWig", selection=BigWigSelection(gr))
      o <- findOverlaps(gr,x)
      tmp <- max(splitAsList(x$score[to(o)], from(o)))
      if(length(tmp)!=length(gr)){
        y <- rep(0,length(gr))
        y[as.integer(names(tmp))] <- as.numeric(tmp)
      }else{
        y <- as.numeric(tmp)
      }
    }else{
      stop("Unsupported file format ", x)
    }
    y
  })
  toExclude <- unique(unlist(lapply(covs, FUN=function(x){
    which(x>quantile(x,maxCovQuant))
  })))
  if(length(toExclude)>0){
    covs <- lapply(covs, FUN=function(x) x[-toExclude])
  }
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
    data.frame(rank=rank(x,ties.method="random"), fraction.of.highest=x/max(x))
  }), .id="file")
  d1$file <- as.factor(d1$file)
  d2$file <- as.factor(d2$file)
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
#' 
#' @examples 
#' # we use an example bigwig file
#' bwf <- system.file("extdata/example_atac.bw", package="epiwraps")
#' # because most of the file is empty, we'll exclude some of the ranges
#' cs <- getCovStats(bwf, exclude=GRanges("1", IRanges(1, 4300000)))
#' plotCovStats(cs)
plotCovStats <- function(qc, labels="AUTO", show.legend=TRUE){
  p1 <- ggplot(qc$coverage, aes(coverage, fraction.above, colour=file)) + 
    geom_line() + labs(x="Read density", y="Fraction of regions > density") +
    theme(legend.position="top", legend.direction = "horizontal") +
    xlim(0,min(max(qc$coverage$coverage),50))
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
  o <- tryCatch( as.dendrogram(hclust(dist(do.call(cbind, dat)))),
                 error=function(e) FALSE )
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


#' fragSizesDist
#'
#' @param x A (named) vector of paths to bam files.
#' @param what Either a positive integer (length 1) indicating how many reads 
#' to randomly sample, or a character vector (of length 1) indicating which
#' chromosome to read.
#' @param flags A `scanBamFlag` object (see \link[Rsamtools]{ScanBamParam})
#' @param BPPARAM A \link[BiocParallel]{BiocParallel} BPPARAM object for
#' multithreading.
#' @param returnPlot Logical; whether to return a plot.
#'
#' @return If `returnPlot=TRUE`, returns a ggplot object, otherwise a
#' data.frame of fragment lengths.
#' @export
#' @importFrom Rsamtools BamFile ScanBamParam scanBamFlag
#' @importFrom GenomicFiles reduceByYield REDUCEsampler
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam SnowParam
#' @importFrom ggplot2 ggplot aes geom_density xlab
#' @import GenomicRanges
#'
#' @examples
#' # example bam file:
#' bam <- system.file("extdata", "ex1.bam", package="Rsamtools")
#' fragSizesDist(bam, what=100)
fragSizesDist <- function(x, what=10000, flags=scanBamFlag(isProperPair=TRUE),
                          BPPARAM=SerialParam(), returnPlot=TRUE){
  stopifnot(length(what)==1)
  if(is.null(names(x))) names(x) <- .cleanFileNames(x)
  flen <- bplapply(x, BPPARAM=BPPARAM, FUN=function(x){
    if(is.numeric(what)){
      what <- as.integer(what)
      stopifnot(what>1)
      x <- GenomicFiles::reduceByYield(
        BamFile(x), MAP=identity, YIELD=function(x)
           readGAlignmentPairs(x, param = ScanBamParam(flag=flags)),
        REDUCE=GenomicFiles::REDUCEsampler(what, verbose=FALSE))
    }else{
      x <- readGAlignmentPairs(x, param=ScanBamParam(flag=flags, 
                               which=GRanges(what, IRanges(1L,10^8))))
    }
    width(as(x, "GRanges"))
  })
  d <- data.frame(sample=rep(factor(names(flen)), lengths(flen)),
                  length=unlist(flen), row.names=NULL)
  if(!returnPlot) return(d)
  ggplot(d, aes(length, colour=sample)) + geom_density() + 
    xlab("Fragment length")
}