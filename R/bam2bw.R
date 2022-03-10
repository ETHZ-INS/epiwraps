#' bam2bw
#' 
#' Creates a coverage bigwig file from a bam file.
#' Note that the implementation requires loading the reads in memory, which can 
#' be quite memory-hungry. Increasing `splitByChr` will decrease the memory 
#' footprint.
#'
#' @param bamfile The path to the bam file.
#' @param output_bw The path to the output bigwig file
#' @param paired Logical; whether to consider fragments instead of reads for 
#'  paired-end data. For spliced (e.g. RNA-seq) data, paired should always be 
#'  set to FALSE.
#' @param binWidth The window size. A lower value (min 1) means a higher 
#'  resolution, but larger file size.
#' @param extend How much to extend single-end reads
#' @param scaling Either TRUE (performs Count Per Million scaling), FALSE (no 
#'   scaling), or a numeric value by which the signal will be divided.
#' @param type Type of the coverage to compile. Either full (full read/fragment),
#'   start (count read/fragment start locations), end, or center.
#' @param strand Strand(s) to capture
#' @param filter Minimum count for a bin to be considered.
#' @param shift Shift (from 3' to 5') by which reads/fragments will be shifted.
#' @param includeDuplicates Logical, whether to include reads flagged as 
#'   duplicates.
#' @param includeSecondary Logical; whether to include secondary alignments
#' @param minMapq Minimum mapping quality (1 to 255)
#' @param minFragLength Minimum fragment length (ignored if `paired=FALSE`)
#' @param maxFragLength Maximum fragment length (ignored if `paired=FALSE`)
#' @param log1p Whether to log-transform (`log(x+1)`) the (scaled) signal.
#' @param splitByChr Whether to process chromosomes separately, and if so by how
#'   many chunks. The should not affect the output, and is simply slightly 
#'   slower and consumes less memory. Can be a logical value, but we instead
#'   recommend giving a positive integer indicating the number of chunks.
#' @param ... Passed to `ScanBamParam`
#' 
#' @return The bigwig filepath
#' @export
#' 
#' @importFrom Rsamtools scanBamHeader scanBamFlag ScanBamParam
#' @importFrom GenomicAlignments readGAlignmentPairs isProperPair readGAlignments
#' @importFrom GenomicRanges resize countOverlaps width score score<- coverage 
#' @importFrom GenomicRanges tileGenome shift
#' @importFrom S4Vectors metadata metadata<-
#' 
#' @examples 
#' # get an example bam file
#' bam <- system.file("extdata", "ex1.bam", package="Rsamtools")
#' # create bigwig
#' bam2bw(bam, scaling=TRUE)
bam2bw <- function(bamfile, output_bw, paired=NULL, binWidth=20L, extend=0L, 
                   scaling=TRUE, type=c("full","center","start","end"),
                   strand=c("*","+","-"), filter=1L, shift=0L, log1p=FALSE,
                   includeDuplicates=TRUE, includeSecondary=TRUE, minMapq=1L, 
                   minFragLength=1L, maxFragLength=5000L, splitByChr=3, 
                   keepSeqLvls=NULL, ...){
  # check inputs
  stopifnot(length(bamfile)==1 && file.exists(bamfile))
  stopifnot(length(output_bw)==1 && is.character(output_bw))
  strand <- match.arg(strand)
  type <- match.arg(type)
  binWidth <- as.integer(binWidth)
  stopifnot(length(binWidth)==1 && binWidth>=1)
  stopifnot(length(filter)==1 && filter>=0)
  if(is.null(paired)){
    message("`paired` not specified, assuming single-end reads.")
    paired <- FALSE
  }
  if(paired) extend <- 0L
  
  seqs <- Rsamtools::scanBamHeader(bamfile)[[1]]$targets
  if(!is.null(keepSeqLvls)){
    if(length(missingLvls <- setdiff(keepSeqLvls), names(seqs))>0)
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

  # obtain coverage (and total count)
  # need to replace the `supressWarnings` with a fix for out-of-chr warnings due
  #   to extension
  res <- suppressWarnings(lapply(names(param), FUN=function(x){
    bam <- .bam2bwGetReads(bamfile, paired, param=param[[x]], scaling=scaling, 
                           type=type, extend=extend, minFragLength=minFragLength,
                           maxFragLength=maxFragLength)
    co <- .bam2bwGetCov(bam, binWidth=binWidth, shift=shift, seqs=chrGroup[[x]], 
                        filter=filter)
    metadata(co)$reads <- metadata(bam)$reads
    rm(bam)
    gc(verbose=FALSE)
    co
  }))

  if(isTRUE(scaling))
    scaling <- sum(unlist(lapply(res, FUN=function(x) metadata(x)$reads )),
                          na.rm=TRUE)/10^6
  
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
  rtracklayer::export.bw(res, output_bw)
  return(invisible(output_bw))
}



.bam2bwGetReads <- function(bamfile, paired, param, scaling, type, extend,
                            minFragLength, maxFragLength){
  if(paired){
    bam <- readGAlignmentPairs(bamfile, param=param)
    bam <- as(bam[isProperPair(bam)], "GRanges")
    if(isTRUE(scaling)) scaling <- length(bam)/10^6
    fl <- width(bam)
    ls <- length(bam)
    bam <- bam[fl>=minFragLength & fl<=maxFragLength]
  }else{
    bam <- as(readGAlignments(bamfile, param=param), "GRanges")
    ls <- length(bam)
    if(type %in% c("full","center") && extend!=0L)
      bam <- resize(bam, width(bam)+as.integer(extend), use.names=FALSE)
  }
  if(type!="full") bam <- resize(bam, width=1L, fix=type, use.names=FALSE)
  metadata(bam)$reads <- ls
  bam
}

.bam2bwGetCov <- function(bam, binWidth, shift, seqs, filter){
  if(binWidth==1L) return(coverage(bam, shift=shift)[names(seqs)])
  if(shift!=0) bam <- shift(bam, shift=shift, use.names=FALSE)
  tiles <- tileGenome(seqs, tilewidth=binWidth, cut.last.tile.in.chrom=TRUE)
  score(tiles) <- countOverlaps(tiles, bam, minoverlap=1L)
  tiles <- tiles[score(tiles)>=as.integer(filter)]
  tiles
}