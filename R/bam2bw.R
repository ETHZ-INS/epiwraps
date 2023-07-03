#' bam2bw
#' 
#' Creates a coverage bigwig file from an alignment bam file.
#'
#' @param bamfile The path to the signal bam file.
#' @param output_bw The path to the output bigwig file
#' @param bgbam The optional path to the control/input bam file to compute 
#'   relative signal (see the `compType` argument).
#' @param paired Logical; whether to consider fragments instead of reads for 
#'  paired-end data. For spliced (e.g. RNA-seq) data, paired should always be 
#'  set to FALSE.
#' @param binWidth The window size. A lower value (min 1) means a higher 
#'  resolution, but larger file size.
#' @param extend The amount *by* which to extend single-end reads (e.g. 
#'   fragment length minus read length). Ignored if `paired=TRUE`.
#' @param scaling Either TRUE (performs Count Per Million scaling), FALSE (no 
#'   scaling), or a numeric value by which the signal will be divided. If 
#'   `bgbam` is given and `scaling=TRUE`, the background will be scaled to the
#'   main signal.
#' @param type Type of the coverage to compile. Either full (full read/fragment),
#'   start (count read/fragment start locations), end, center, or 'ends' (both
#'   ends of the read/fragment).
#' @param strand Strand(s) to capture (any by default).
#' @param strandMode The strandMode of the data (whether the strand is given by
#'   the first or second mate, which depends on the library prep protocol). See
#'   \link[GenomicAlignments]{strandMode} for more information. This parameter 
#'   has no effect unless one of the `strand`, `extend` parameters or a 
#'   strand-specific `shift` are used.
#' @param shift Shift (from 3' to 5') by which reads/fragments will be shifted.
#'   If `shift` is an integer vector of length 2, the first value will represent
#'   the shift for the positive strand, and the second for the negative strand.
#' @param includeDuplicates Logical, whether to include reads flagged as 
#'   duplicates.
#' @param includeSecondary Logical; whether to include secondary alignments
#' @param minMapq Minimum mapping quality (1 to 255)
#' @param minFragLength Minimum fragment length (ignored if `paired=FALSE`)
#' @param maxFragLength Maximum fragment length (ignored if `paired=FALSE`)
#' @param log1p Whether to log-transform (`log(x+1)`) the (scaled) signal.
#'   Ignored when the signal is relative to a background.
#' @param keepSeqLvls An optional vector of seqLevels (i.e. chromosomes) to 
#'   include.
#' @param splitByChr Whether to process chromosomes separately, and if so by how
#'   many chunks. The should not affect the output, and is simply slightly 
#'   slower and consumes less memory. Can be a logical value (in which case each
#'   chromosome is processed separately), but we instead recommend giving a 
#'   positive integer indicating the number of chunks.
#' @param compType The type of relative signal to compute (ignored if `bgbam` 
#'   isn't given). Can either be 'log2ratio' (log2-foldchange between signals),
#'   'subtract' (difference), 'log10ppois' (rounded -log10 poisson p-value),
#'   'log10FE' (log10 fold-enrichment). For fold-based signals, `pseudocount` 
#'   will be added before computing the ratio. By default negative values are
#'   capped (see `zeroCap`).
#' @param zeroCap Logical; whether to cap values below zero to zero when 
#'   producing relative signals.
#' @param pseudocount The count to be added before fold-enrichment calculation 
#'   between signals. Ignored if `compType="subtract"`.
#' @param localBackground A (vector of) number of windows around each 
#'   tile/position for which a local background will be calculated. The 
#'   background used will be the maximum of the one at the given position or of
#'   the mean of each of the local backgrounds.
#' @param forceSeqlevelsStyle If specified, forces the use of the specified 
#'   seqlevel style for the output bigwig. Can take any value accepted by
#'   `seqlevelsStyle`.
#' @param exclude An optional GRanges of regions for which overlapping reads 
#'   should be excluded.
#' @param only An optional GRanges of regions for which overlapping reads should
#'   be included. If set, all other reads are discarded.
#' @param binSummarization The method to summarize nucleotides into each bin,
#'   either "max" (default), "min" or "mean".
#' @param verbose Logical; whether to print progress messages
#' @param ... Passed to `ScanBamParam`
#' 
#' @return The bigwig filepath. Alternatively, if `output_bw=NA_character_`, 
#'   the coverage data is not written to file but returned.
#' 
#' @details 
#' * For single-end ChIPseq data, extend reads to the expected fragment size 
#'   using the `extend` argument.
#' * The implementation involves reading the reads before processing, which can
#'  be quite memory-hungry. To avoid this the files are read by chunks (composed
#'   of chromosomes of balanced sizes), controlled by the `splitByChr` argument.
#'   This reduces memory usage but also creates overhead which reduces the 
#'   speed. In contexts where the file is small or memory isn't a problem, using 
#'  `splitByChr=FALSE` will achieve the best speed. Otherwise, increasing 
#'  `splitByChr` will decrease the memory footprint.
#' * Consider restricting the read quality using `includeDuplicates=FALSE` 
#'   and `minMapq=20` for example.
#' 
#' @author Pierre-Luc Germain
#' 
#' @export
#' 
#' @importFrom Rsamtools scanBamHeader scanBamFlag ScanBamParam
#' @importFrom GenomicAlignments readGAlignmentPairs isProperPair readGAlignments
#' @importFrom GenomicRanges resize countOverlaps width score score<- coverage 
#' @importFrom GenomicRanges tileGenome shift trim
#' @importFrom GenomeInfoDb Seqinfo seqinfo seqinfo<-
#' @importFrom S4Vectors metadata metadata<- runmean Rle
#' @importFrom IRanges RleList
#' @importFrom pbapply pblapply
#' 
#' @examples 
#' # get an example bam file
#' bam <- system.file("extdata", "ex1.bam", package="Rsamtools")
#' # create bigwig
#' bam2bw(bam, "output.bw")
bam2bw <- function(bamfile, output_bw, bgbam=NULL, paired=NULL,
                   binWidth=20L, extend=0L, scaling=TRUE, shift=0L,
                   type=c("full","center","start","end","ends"),
                   compType=c("log2ratio", "subtract", "log10ppois", "log10FE"),
                   strand=c("*","+","-"), strandMode=1, log1p=FALSE, exclude=NULL,
                   includeDuplicates=TRUE, includeSecondary=FALSE, minMapq=1L, 
                   minFragLength=1L, maxFragLength=5000L, keepSeqLvls=NULL, 
                   splitByChr=3, pseudocount=1L, localBackground=1L, only=NULL,
                   zeroCap=TRUE, forceSeqlevelsStyle=NULL, verbose=TRUE, 
                   binSummarization=c("max","min","mean"), ...){
  # check inputs
  stopifnot(length(bamfile)==1 && file.exists(bamfile))
  if(!is.null(exclude)) stopifnot(is(exclude,"GRanges"))
  if(!is.null(bgbam)){
    stopifnot(length(bgbam)==1 && file.exists(bgbam))
    stopifnot(is.logical(scaling) && length(scaling)==1)
  }
  stopifnot(length(output_bw)==1 &&
              (is.character(output_bw) || is.na(output_bw)))
  compType <- match.arg(compType)
  strand <- match.arg(strand)
  type <- match.arg(type)
  binWidth <- as.integer(binWidth)
  shift <- as.integer(shift)
  stopifnot(length(shift) %in% 1:2)
  stopifnot(length(binWidth)==1 && binWidth>=1)
  if(is.null(paired)){
    if(verbose) message("`paired` not specified, assuming single-end reads. ",
                        "Set to paired='auto' to automatically detect.")
    paired <- FALSE
  }else if(paired=="auto"){
    paired <- testPairedEndBam(bamfile)
    if(verbose) message("Detected ", ifelse(paired,"paired","unpaired")," data")
  }
  if(paired) extend <- 0L
  if(type=="ends" && !paired && verbose)
    warning("type='ends' only makes sense with paired-end libraries...")

  # prepare flags for bam reading
  strandflg <- switch(strand, "*"=NA, "+"=FALSE, "-"=TRUE)
  if(paired) strandflg <- NA
  flgs <- scanBamFlag(isDuplicate=ifelse(includeDuplicates,NA,FALSE), 
                      isSecondaryAlignment=ifelse(includeSecondary,NA,FALSE),
                      isMinusStrand=strandflg,isNotPassingQualityControls=FALSE)  
    
  seqs <- Rsamtools::scanBamHeader(bamfile)[[1]]$targets
  param <- .getBamChunkParams(bamfile, flgs=flgs, keepSeqLvls=keepSeqLvls, 
                              nChunks=splitByChr)
  seqs <- seqs[.checkMissingSeqLevels(seqs, keepSeqLvls)]
  
  if(verbose) message("Reading in signal...")
  
  # obtain coverage (and total count)
  res <- pblapply(names(param), FUN=function(x){
    .bam2bwReadChunk(bamfile, param=param[[x]], binWidth=binWidth,
                     paired=paired, type=type, extend=extend, shift=shift, 
                     minFragL=minFragLength, maxFragL=maxFragLength, 
                     forceStyle=forceSeqlevelsStyle, exclude=exclude,
                     keepStrand=ifelse(paired && strand!="*",strand,"*"),
                     binSummarization=binSummarization, si=seqs,
                     strandMode=strandMode, only=only)
  })
  
  if(is.null(bgbam)){
    if(verbose) message("Writing bigwig...")
    res <- .bam2bwScaleCovList(res, scaling=scaling, log1p=log1p)
    if(is.na(output_bw)) return(res)
    rtracklayer::export.bw(res, output_bw)
    return(invisible(output_bw))
  }
    
  # from here on: relative signal
  
  if(verbose) message("Reading in background...")
  
  # same as above
  bg <- pblapply(names(param), FUN=function(x){
    .bam2bwReadChunk(bgbam, param=param[[x]], binWidth=binWidth,
                     paired=paired, type=type, extend=extend, shift=shift, 
                     minFragL=minFragLength, maxFragL=maxFragLength, 
                     forceStyle=forceSeqlevelsStyle, exclude=exclude,
                     keepStrand=ifelse(paired && strand!="*",strand,"*"),
                     binSummarization=binSummarization, si=seqs,
                     strandMode=strandMode, only=only)
  })

  if(verbose) message("Computing relative signal...")
  
  res <- .bam2bwScaleCovList(res, scaling=FALSE)
  bg <- .bam2bwScaleCovList(bg, scaling=FALSE)
  if(isTRUE(scaling)){
    bg <- bg * .covTrimmedMean(res)/.covTrimmedMean(bg)
  }
  
  res <- .bwRelativeSignal(res, bg, compType=compType, pseudocount=pseudocount,
                           localBackground=localBackground, zeroCap=zeroCap)
  rm(bg)
  
  if(is.na(output_bw)) return(res)
  
  if(verbose) message("Writing bigwig...")
  rtracklayer::export.bw(res, output_bw)
  return(invisible(output_bw))
}

# expects a single Rle, not RleList
.bwRelativeSignal <- function(res, bg, pseudocount=1L, localBackground=1L,
                              compType=c("log2ratio", "subtract", 
                                         "log10ppois", "log10FE"),
                              zeroCap=TRUE){
  compType <- match.arg(compType)
  res <- lapply(setNames(names(res), names(res)), FUN=function(x){
    bg <- .bam2bwLocalBackground(bg[[x]], windows=localBackground)
    res <- res[[x]]
    if(compType=="log10ppois"){
      res[res<=pseudocount] <- 0L
      y <- Rle(as.integer(round(-log10(ppois(
            as.integer(res), as.numeric(bg)+pseudocount, lower.tail=FALSE)))))
    }else{
      y <- switch(compType,
                  log2ratio=log10((res+pseudocount)/(bg+pseudocount)),
                  subtract=res-bg,
                  logFE=log10((res+pseudocount)/(bg+pseudocount)),
            )
      if(zeroCap) y[y<0L] <- 0L
    }
    y
  })
  as(res, "RleList")
}


#' frag2bw
#' 
#' Creates a coverage bigwig file from a Tabix-indexed fragment file.
#'
#' @param tabixFile The path to a tabix-indexed bam file, or a TabixFile object.
#' @param output_bw The path to the output bigwig file
#' @param barcodes An optional list of barcodes to use (assuming that the file
#'   contains the column)
#' @param binWidth The window size. A lower value (min 1) means a higher 
#'  resolution, but larger file size.
#' @param scaling Either TRUE (performs Count Per Million scaling), FALSE (no 
#'   scaling), or a numeric value by which the signal will be divided. If 
#'   `bgbam` is given and `scaling=TRUE`, the background will be scaled to the
#'   main signal.
#' @param type Type of the coverage to compile. Either full (full read/fragment),
#'   start (count read/fragment start locations), end, center, or 'ends' (both
#'   ends of the read/fragment).
#' @param strand Strand(s) to capture (any by default).
#' @param shift Shift (from 3' to 5') by which reads/fragments will be shifted.
#'   If `shift` is an integer vector of length 2, the first value will represent
#'   the shift for the positive strand, and the second for the negative strand.
#' @param useScore Whether to use the score column (if any) as coverage weights.
#' @param minFragLength Minimum fragment length (ignored if `paired=FALSE`)
#' @param maxFragLength Maximum fragment length (ignored if `paired=FALSE`)
#' @param keepSeqLvls An optional vector of seqLevels (i.e. chromosomes) to 
#'   include.
#' @param forceSeqlevelsStyle If specified, forces the use of the specified 
#'   seqlevel style for the output bigwig. Can take any value accepted by
#'   `seqlevelsStyle`.
#' @param exclude An optional GRanges of regions for which overlapping reads 
#'   should be excluded.
#' @param binSummarization The method to summarize nucleotides into each bin,
#'   either "max" (default), "min" or "mean".
#' @param verbose Logical; whether to print progress messages
#' 
#' @return The bigwig filepath. Alternatively, if `output_bw=NA_character_`, 
#'   the coverage data is not written to file but returned.
#' 
#' @export
frag2bw <- function(tabixFile, output_bw, binWidth=20L, extend=0L, scaling=TRUE,
                    type=c("full","center","start","end","ends"), barcodes=NULL,
                    strand=c("*","+","-"), shift=0L, log1p=FALSE, exclude=NULL,
                    minFragLength=1L, maxFragLength=5000L, keepSeqLvls=NULL, 
                    useScore=FALSE, forceSeqlevelsStyle=NULL, 
                    binSummarization=c("max","min","mean"), verbose=TRUE, ...){
  binSummarization <- match.arg(binSummarization)
  if(!is(tabixFile, "TabixFile")) tabixFile <- TabixFile(tabixFile)
  
  if(verbose) message("Reading in signal...")
  res <- tabixChrApply(tabixFile, keepSeqLvls=keepSeqLvls, exclude=exclude,
                       FUN=function(x){
    if(!is.null(barcodes) && !is.null(x$name))
      x <- x[which(x$name %in% barcodes)]
    .frag2co(x, strand=strand, type=type, minFragLength=minFragLength,
             maxFragLength=maxFragLength, shift=shift, useScore=useScore,
             binSummarization=binSummarization)
  })
  if(verbose) message("Writing bigwig...")
  res <- .bam2bwScaleCovList(res, scaling=scaling)
  if(is.na(output_bw)) return(res)
  rtracklayer::export.bw(res, output_bw)
  return(invisible(output_bw))
}

.frag2co <- function(x, strand, minFragLength, maxFragLength, type, shift,
                     binSummarization="max", useScore=FALSE){
  nr <- length(x)
  if(strand!="*") x <- x[strand(x)==strand]
  x <- x[width(x)>=minFragLength & width(x)<=maxFragLength]
  x <- .doShift(x, shift)
  if(type=="ends"){
    x <- .align2cuts(x)
  }else if(type!="full"){
    x <- resize(x, width=1L, fix=type, use.names=FALSE)
  }
  x <- coverage(x, weight=ifelse(useScore && !is.null(x$score),"score",1L))
  co <- tileRle(x, bs=binWidth, method=binSummarization)
  metadata(co)$reads <- nr
  rm(x)
  gc(full=TRUE, verbose=FALSE)
  co
}

#' @importFrom GenomeInfoDb seqlevelsStyle<- seqlevelsInUse
.bam2bwReadChunk <- function(bamfile, param, binWidth, forceStyle=NULL, si=NULL,
                             keepStrand="*", binSummarization="max", only=NULL, 
                             ...){
  # get reads/fragments from chunk
  bam <- .bam2bwGetReads(bamfile, param=param, si=si, only=only,
                         forceStyle=forceStyle, ...)
  if(keepStrand != "*") bam <- bam[which(strand(bam)==keepStrand)]
  # compute coverages
  co <- tileRle(coverage(bam), bs=binWidth, method=binSummarization)
  # save library size for later normalization
  metadata(co)$reads <- metadata(bam)$reads
  rm(bam)
  gc(full=TRUE, verbose=FALSE)
  co
}


# reads reads from bam file.
.bam2bwGetReads <- function(bamfile, paired, param, type, extend, shift=0L,
                            minFragL, maxFragL, forceStyle=NULL, si=NULL,
                            only=NULL, exclude=NULL, strandMode=0){
  if(paired){
    bam <- readGAlignmentPairs(bamfile, param=param, strandMode=strandMode)
    bam <- as(bam[isProperPair(bam)], "GRanges")
    ls <- length(bam)
    bam <- bam[width(bam)>=minFragL & width(bam)<=maxFragL]
  }else{
    bam <- as(readGAlignments(bamfile, param=param), "GRanges")
    ls <- length(bam)
    if(type %in% c("full","center") && extend!=0L)
      bam <- suppressWarnings(trim(resize(bam, width(bam)+as.integer(extend),
                                          use.names=FALSE)))
  }
  if(!is.null(forceStyle)) seqlevelsStyle(bam) <- forceStyle
  if(!is.null(only) && is(only,"GRanges")){
    .comparableStyles(bam, only)
    bam <- bam[overlapsAny(bam,only)]
  }
  if(!is.null(exclude) && is(exclude,"GRanges")){
    .comparableStyles(bam, exclude)
    bam <- bam[!overlapsAny(bam,exclude)]
  }
  bam <- .doShift(bam, shift)
  if(type=="ends"){
    bam <- .align2cuts(bam)
  }else if(type!="full"){
    bam <- resize(bam, width=1L, fix=type, use.names=FALSE)
  }
  if(!is.null(si)){
    # for trimming, to avoid out-of-range warnings
    bam <- keepSeqlevels(bam, names(si), pruning.mode="coarse")
    if(is.vector(si)) si <- Seqinfo(names(si), as.integer(si))
    seqinfo(bam) <- si
  }
  bam <- trim(bam)
  bam <- bam[width(bam)>0]
  metadata(bam)$reads <- ls
  bam
}

.doShift <- function(x, shift){
  stopifnot(is(x, "GRanges"))
  shift <- as.integer(shift)
  stopifnot(length(shift) %in% 1:2)
  if(all(shift==0)) return(x)
  if(length(shift)>1){
    w <- strand(x)=="+"
    bam <- sort(c(
      shift(x[IRanges::which(w)], shift=shift[1], use.names=FALSE),
      shift(x[IRanges::which(!w)], shift=shift[2], use.names=FALSE)))
  }else{
    x <- shift(x, shift=shift, use.names=FALSE)
  }
  x
}

#' tileRle
#' 
#' Creates an Rle of fixed-with bins from a continuous numeric Rle
#'
#' @param x A numeric `Rle` (or `RleList`)
#' @param bs A positive integer specifying the bin size
#' @param method The method for summarizing bins
#' @param roundSummary Logical; whether to round bins with summarized coverage
#'   (default FALSE)
#'
#' @return An object of same class and length as `x`
#' @export
#' @importFrom IRanges viewMins median quantile
#'
#' @examples
#' # creating a dummy coverage and visualizing it:
#' cov <- Rle(rpois(100,0.5))
#' plot(cov, type="l", col="lightgrey")
#' # summarizing to tiles of width 5 (by default using maximum)
#' cov2 <- tileRle(cov, bs=5L)
#' lines(cov2, col="red")
tileRle <- function(x, bs=10L, method=c("max","min","mean"), roundSummary=FALSE){
  bs <- as.integer(bs)
  stopifnot(bs>=1L)
  min(min(runLength(x)[lengths(runLength(x))>0]))
  if(bs<=min(min(runLength(x)[lengths(runLength(x))>0]))) return(x)
  method <- match.arg(method)
  if(is(x,"RleList"))
    return(as(lapply(x, bs=bs, method=method, FUN=tileRle), "RleList"))
  stopifnot(is(x,"Rle"))
  
  # define the new ends of runs based on the number of full 'bs'
  cs <- cumsum(runLength(x)) # gives the end position of each run
  # floored end tiled position, except for the last bit (to respect chr sizes)
  csf <- c(bs*floor(cs[-length(cs)]/bs), cs[length(cs)])
  # new run lengths based on the difference to the floored end of previous run
  new.lengths <- pmax(c(csf[-length(csf)],tail(cs,1)) - 
                        (bs*ceiling(c(0L,cs[-length(cs)])/bs)), 0L)

  # identify the runs whose ends doesn't fall on a tile end
  w <- which(cs > csf)
  # get the values for tiles that overlap more than one run
  bins <- Views(x, IRanges(unique(csf[w])+1L, width=bs))
  binVals <- switch(method,
                    max=viewMaxs(bins),
                    min=viewMins(bins),
                    mean=viewMeans(bins))
  if(roundSummary) binVals <- round(binVals)
  
  # remove runs that are entirely contained in a tile (except very last)
  removedRuns <- setdiff(which(runLength(x)<bs), length(cs))
  new.lengths[removedRuns] <- 0L
  
  # create new runs for boundary tiles, and inject them in the Rle vectors
  extra.w <- w[!duplicated(csf[w])]
  Rle(values=inject(binVals, runValue(x), at=extra.w),
      lengths=inject(bs, pmax(new.lengths,0L), at=extra.w))
}


# flattens a list of coverages (RleList) or of scored windows (GRangesList) and 
# eventually applies a scaling
.bam2bwScaleCovList <- function(res, scaling, log1p=FALSE){
  if(isTRUE(scaling)){
    # get scaling from chunks read counts
    scaling <- sum(unlist(lapply(res, FUN=function(x) metadata(x)$reads )),
                   na.rm=TRUE)/10^6
    if(is.na(scaling) || !isTRUE(scaling>0))
      stop(paste("Cannot establish scaling factor; object most likely doesn't",
                 "contain a metadata(x)$reads."))
  }
  if(is(res[[1]], "GRanges")){
    res <- unlist(GRangesList(res))
  }else{
    res <- .reduceRleLists(res)
  }
  if(!isFALSE(scaling)){
    stopifnot(scaling>0)
    if(is(res, "RleList")){
      res <- res/scaling
    }else{
      score(res) <- score(res)/scaling
    }
  }
  if(log1p){
    if(is(res, "RleList")){
      res <- log1p(res)
    }else{
      score(res) <- log1p(score(res))
    }
  }
  res
}

.reduceRleLists <- function(res, fn="+"){
  sn <- unique(unlist(lapply(res, names)))
  as(lapply(setNames(sn,sn), FUN=function(x){
    r2 <- lapply(res, FUN=function(co){
      if(x %in% names(co)) return(co[[x]])
      return(NULL)
    })
    suppressWarnings(Reduce(fn, r2[which(!sapply(r2,is.null))]))
  }), "RleList")
}

# calculates local background at positions, analogous to MACS
# TO DO: handle chr that are smaller than windows !!!
.bam2bwLocalBackground <- function(x, windows=10L){
  if(length(windows)==0 || sum(windows)<=1L) return(x)
  if(is(x,"GRanges")){
    x <- x[order(seqnames(x))]
    score(x) <- Rle(score(x))
    score(x) <- unlist(.bam2bwLocalBackground(
      split(Rle(score(x)), seqnames(x)), windows=windows))
    return(x)
  }
  if(is.numeric(x)) x <- Rle(x)
  stopifnot(is(x,"RleList") || is(x,"Rle"))
  if(length(windows)==1) return(runmean(x, k=windows, endrule = "constant"))
  do.call(pmax, lapply(windows, FUN=function(w){
    if(w==1) return(x)
    runmean(x, k=w, endrule = "constant")
  }))
}

.align2cuts <- function(x){
  x <- granges(GRanges(x))
  sort(c(resize(x, width=1L, fix="start", use.names=FALSE),
         resize(x, width=1L, fix="end", use.names=FALSE)))
}

# returns the bin size of a fixed-width RleList
.getBinWidth <- function(x){
  if(is(x,"RleList")){
    x <- unlist(lapply(x, .getBinWidth))
    return(as.integer(round(median(x, na.rm=TRUE))))
  }else if(is(x,"Rle")){
    x <- runLength(x)
    if(length(x)==1) return(NA_integer_)
    return(min(x[-length(x)]))
  }
  stop("x is not a Rle/RleList")
}
