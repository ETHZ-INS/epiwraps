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
#' @param filter Minimum count (in the signal) for a bin to be considered.
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
#'   isn't given). Can either be `ppois` (-log10 Poisson p-value), `logFE` 
#'   (log10 fold-enrichment), or `subtract`. If `logFE`, the `pseudocount` 
#'   will be added before computing the ratio.
#' @param pseudocount The count to be added before fold-enrichment calculation 
#'   between signals. Only used for `compType="logFE"`.
#' @param localBackground Whether to use the local background instead of the 
#'   count at the specific position/window when computing relative signal 
#'   (analogous to MACS). Can either be a logical value, or an integer vector of
#'   multiple sizes at which to compute a local background (as average signal).
#'   The maximum of the average background signal across the window sizes will 
#'   be used. If TRUE, will use 1kb and 5kb.
#' @param verbose Logical; whether to print progress messages
#' @param ... Passed to `ScanBamParam`
#' 
#' @return The bigwig filepath
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
#' 
#' @examples 
#' # get an example bam file
#' bam <- system.file("extdata", "ex1.bam", package="Rsamtools")
#' # create bigwig
#' bam2bw(bam, "output.bw")
bam2bw <- function(bamfile, output_bw, bgbam=NULL, paired=NULL, binWidth=20L, 
                   extend=0L, compType=c("ppois", "subtract", "logFE"),
                   scaling=TRUE, type=c("full","center","start","end","ends"),
                   strand=c("*","+","-"), filter=1L, shift=0L, log1p=FALSE,
                   includeDuplicates=TRUE, includeSecondary=FALSE, minMapq=1L, 
                   minFragLength=1L, maxFragLength=5000L, keepSeqLvls=NULL, 
                   splitByChr=3, pseudocount=1L, localBackground=c(1000L,5000L),
                   verbose=TRUE, ...){
  # check inputs
  stopifnot(length(bamfile)==1 && file.exists(bamfile))
  if(!is.null(bgbam)){
    stopifnot(length(bgbam)==1 && file.exists(bgbam))
    stopifnot(is.logical(scaling) && length(scaling)==1)
  }
  stopifnot(length(output_bw)==1 && is.character(output_bw))
  compType <- match.arg(compType)
  strand <- match.arg(strand)
  type <- match.arg(type)
  binWidth <- as.integer(binWidth)
  shift <- as.integer(shift)
  stopifnot(length(shift) %in% 1:2)
  stopifnot(length(binWidth)==1 && binWidth>=1)
  stopifnot(length(filter)==1 && filter>=0)
  if(is.null(paired)){
    if(verbose) message("`paired` not specified, assuming single-end reads.")
    paired <- FALSE
  }
  if(paired) extend <- 0L
  if(type=="ends" && !paired && verbose)
    warning("type='ends' only makes sense with paired-end libraries...")
  
  seqs <- Rsamtools::scanBamHeader(bamfile)[[1]]$targets
  if(!is.null(keepSeqLvls)){
    if(length(missingLvls <- setdiff(keepSeqLvls, names(seqs))>0))
      stop(paste0(
        "Some of the seqLevels specified by `keepSeqLvls` are not in the data.
The first few are:",
        head(paste(missingLvls, collapse=", "), 3)))
    seqs <- seqs[keepSeqLvls]
  }
  
  # prepare flags for bam reading
  flgs <- scanBamFlag(isDuplicate=ifelse(includeDuplicates,NA,FALSE), 
                      isSecondaryAlignment=ifelse(includeSecondary,NA,FALSE),
                      isMinusStrand=switch(strand, "*"=NA, "+"=FALSE, "-"=TRUE),
                      isNotPassingQualityControls=FALSE)
  
  # generate list of reading params
  if(splitByChr){
    if(is.numeric(splitByChr) && splitByChr>=2){
      splitByChr <- as.integer(splitByChr)
      seqs <- sort(seqs, dec=TRUE)
      chrGroup <- head(rep(seq_len(splitByChr), 
                           ceiling(length(seqs)/splitByChr)), length(seqs))
      chrGroup <- split(seqs, chrGroup)
    }else{
      chrGroup <- split(seqs, names(seqs))
    }
    param <- lapply(chrGroup, FUN=function(x)
      ScanBamParam(flag=flgs, mapqFilter=minMapq, 
                   which=GRanges(names(x), IRanges(1L,x)), ...))
  }else{
    chrGroup <- list(wg=seqs)
    param <- list(wg=ScanBamParam(flag=flgs, mapqFilter=minMapq, ..., 
                                  which=GRanges(names(seqs), IRanges(1L,seqs))))
  }

  if(verbose) message("Reading in signal...")
  
  # obtain coverage (and total count)
  # need to replace the `supressWarnings` with a fix for out-of-chr warnings due
  #   to read extension
  res <- lapply(names(param), FUN=function(x){
    .bam2bwReadChunk(bamfile, param=param[[x]], binWidth=binWidth,
                     paired=paired, type=type, extend=extend, shift=shift, 
                     minFragL=minFragLength, maxFragL=maxFragLength, 
                     seqs=chrGroup[[x]])
  })
  
  if(is.null(bgbam)){
    if(verbose) message("Writing bigwig...")
    res <- .bam2bwScaleCovList(res, scaling=scaling, log1p=log1p)
    rtracklayer::export.bw(res, output_bw)
    return(invisible(output_bw))
  }
    
  # from here on: relative signal
  
  if(verbose) message("Reading in background...")
  
  # same as above, except that we the same tiles and no filter if bin-based
  bg <- lapply(names(param), FUN=function(x){
    tiles <- NULL
    if(is(res[[x]], "GRanges")) tiles <- res[[x]]
    .bam2bwReadChunk(bgbam, param=param[[x]], binWidth=binWidth,
                     paired=paired, type=type, extend=extend, shift=shift, 
                     minFragL=minFragLength, maxFragL=maxFragLength, 
                     seqs=chrGroup[[x]], filter=0, tiles=tiles)
  })

  if(verbose) message("Computing relative signal...")
  
  if(scaling){
    scaling <- sum(unlist(lapply(res, FUN=function(x) metadata(x)$reads )),
                   na.rm=TRUE)/
                sum(unlist(lapply(bg, FUN=function(x) metadata(x)$reads )),
                   na.rm=TRUE)
  }
  res <- .bam2bwScaleCovList(res, scaling=FALSE)
  bg <- .bam2bwScaleCovList(bg, scaling=scaling)
  
  if(!isFALSE(localBackground)){
    if(isTRUE(localBackground)) localBackground <- 1000L*c(1L,5L)
    bg <- .bam2bwLocalBackground(bg, binWidth=binWidth, windows=localBackground)
  }
  
  if(is(res, "RleList")){
    res <- switch(compType,
            subtract=RleList(lapply(setNames(names(res), names(res)), 
                                    FUN=function(x) pmax(res[[x]]-bg[[x]],0))),
            logFE=log10((res+pseudocount)/(bg+pseudocount)),
            ppois=RleList(lapply(setNames(names(res), names(res)), 
                                 FUN=function(x){
              Rle(-log10(ppois(as.integer(res[[x]]), as.numeric(bg[[x]]), 
                        lower.tail=FALSE)))
            }))
          )
  }else{
    score(res) <- switch(compType,
             subtract=pmax(score(res)/score(bg), 0),
             logFE=log10((score(res)+pseudocount)/(score(bg)+pseudocount)),
             ppois=-log10(ppois(score(res), score(bg), lower.tail=FALSE)))
  }
  rm(bg)
  
  if(verbose) message("Writing bigwig...")

  rtracklayer::export.bw(res, output_bw)
  return(invisible(output_bw))
}

.bam2bwReadChunk <- function(bamfile, param, binWidth, seqs, filter=0,
                             tiles=NULL, ...){
  # get reads/fragments from chunk
  bam <- .bam2bwGetReads(bamfile, param=param, si=seqs, ...)
  # compute coverages
  co <- .bam2bwGetCov(bam, binWidth=binWidth, seqs=seqs, 
                      tiles=tiles, filter=filter)
  # save library size for later normalization
  metadata(co)$reads <- metadata(bam)$reads
  rm(bam)
  gc(verbose=FALSE)
  co
}


# reads reads from bam file.
.bam2bwGetReads <- function(bamfile, paired, param, type, extend,
                            shift, minFragL, maxFragL, si=NULL){
  if(paired){
    bam <- readGAlignmentPairs(bamfile, param=param)
    bam <- as(bam[isProperPair(bam)], "GRanges")
    ls <- length(bam)
    bam <- bam[width(bam)>=minFragL & width(bam)<=maxFragL]
  }else{
    bam <- as(readGAlignments(bamfile, param=param), "GRanges")
    ls <- length(bam)
    if(type %in% c("full","center") && extend!=0L)
      bam <- resize(bam, width(bam)+as.integer(extend), use.names=FALSE)
  }
  if(!all(shift==0)){
    if(length(shift)>1){
      bam[strand(bam)=="+"] <- shift(bam[strand(bam)=="+"], shift=shift[1], 
                                     use.names=FALSE)
      bam[strand(bam)=="-"] <- shift(bam[strand(bam)=="-"], shift=shift[2], 
                                     use.names=FALSE)
    }else{
      bam <- shift(bam, shift=shift, use.names=FALSE)
    }
  }
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
  metadata(bam)$reads <- ls
  bam
}

# reads per base or per window coverage from GRanges of fragments
.bam2bwGetCov <- function(bam, binWidth, seqs, filter, tiles=NULL){
  if(binWidth==1L) return(coverage(bam)[names(seqs)])
  if(is.null(tiles))
    tiles <- tileGenome(seqs, tilewidth=binWidth, cut.last.tile.in.chrom=TRUE)
  score(tiles) <- countOverlaps(tiles, bam, minoverlap=1L)
  tiles <- tiles[score(tiles)>=as.integer(filter)]
  tiles
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
    res <- lapply(res, as.list)
    names(res) <- NULL
    res <- do.call(RleList, unlist(res, recursive=FALSE))
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

# calculates local background at positions, analogous to MACS
.bam2bwLocalBackground <- function(x, binWidth=1L, windows=1000L*c(1L,5L)){
  if(is(x,"GRanges")){
    x <- x[order(seqnames(x))]
    score(x) <- Rle(score(x))
    score(x) <- unlist(.bam2bwLocalBackground(
      split(Rle(score(x)), seqnames(x)), binWidth=binWidth, windows=windows))
    return(x)
  }
  if(is.numeric(x)) x <- Rle(x)
  stopifnot(is(x,"RleList") || is(x,"Rle"))
  windows <- windows[windows>=binWidth]
  if(length(windows)==0) return(x)
  windows <- as.integer(windows/binWidth)
  if(length(windows)==1) return(runmean(x, k=windows, endrule = "constant"))
  do.call(pmax, lapply(windows, FUN=function(w){
    if(w==1) return(x)
    runmean(x, k=w, endrule = "constant")
  }))
}
