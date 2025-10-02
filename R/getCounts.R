#' getCounts
#' 
#' Creates a SummarizedExperiment of fragment (or insertion) counts from bam 
#' files that overlap given regions.
#'
#' @param bam_files A vector of paths to the bam files.
#' @param regions A `GRanges` of regions in which to counts.
#' @param paired Logical; whether the data is paired (assumed unpaired by 
#'   default). Use `paired="auto"` for automatic detection using the first bam 
#'   file.
#' @param ignore.strand Logical; whether to ignore strand for the purpose of
#'   counting overlaps (default TRUE).
#' @param randomAcc Logical, whether to use random access. This is disabled by
#'   default because the overhead of random access to a lot of regions is 
#'   typically worse than reading the entire file. However, if you need to get
#'   counts in few regions, enabling this will be faster. Note however that 
#'   when using random access, the output object will not contain depth 
#'   information.
#' @param ov.type Overlap type. See the `type` argument of 
#'   \code{link[GenomicRanges]{countOverlaps}}.
#' @param maxgap Maximum gap allowed for overlaps (see the corresponding 
#'   argument of \code{link[GenomicRanges]{countOverlaps}}).
#' @param minoverlap Minimum overlap (see the corresponding argument of 
#'   \code{link[GenomicRanges]{countOverlaps}}).
#' @inheritParams bam2bw
#' 
#' @return A RangedSummarizedExperiment with a 'counts' assay.
#' @export
#'
#' @examples
#' # get an example bam file
#' bam <- system.file("extdata", "ex1.bam", package="Rsamtools")
#' # create regions of interest
#' peaks <- GRanges(c("seq1","seq1","seq2"), IRanges(c(400,900,500), width=100))
#' getCounts(bam, peaks, paired=FALSE)
getCounts <- function(bam_files, regions, paired, extend=0L, shift=0L, 
                      type=c("full","center","start","end","ends"),
                      ov.type="any", maxgap=-1L, minoverlap=1L,
                      ignore.strand=TRUE, strandMode=1, includeDuplicates=TRUE, 
                      includeSecondary=FALSE, minMapq=1L, minFragLength=1L,
                      maxFragLength=5000L, splitByChr=3, randomAcc=FALSE,
                      verbose=TRUE){
  # check inputs
  stopifnot(is(regions, "GRanges"))
  stopifnot(all(vapply(bam_files,FUN.VALUE=logical(1L),FUN=file.exists)))
  stopifnot(all(grepl("\\.bam",bam_files,ignore.case=TRUE)))
  type <- match.arg(type)
  stopifnot(extend>=0L)
  if(is.null(paired)){
    if(verbose) message("`paired` not specified, assuming single-end reads. ",
                        "Set to paired='auto' to automatically detect.")
    paired <- FALSE
  }else if(paired=="auto"){
    paired <- testPairedEndBam(bamfile)
    if(verbose) message("Detected ", ifelse(paired,"paired","unpaired")," data")
  }
  if(paired && !(type %in% c("center","ends")) && verbose)
    message("`extend` argument ignored.")
  if(type=="ends" && !paired){
    if(verbose)
      warning("type='ends' typically only makes sense with paired-end data...")
    type <- "start"
  }
  # prepare flags for bam reading
  flgs <- scanBamFlag(isDuplicate=ifelse(includeDuplicates,NA,FALSE), 
                      isSecondaryAlignment=ifelse(includeSecondary,NA,FALSE),
                      isNotPassingQualityControls=FALSE)
  seqs <- Rsamtools::scanBamHeader(bam_files[[1]])[[1]]$targets
  seqs <- seqs[.checkMissingSeqLevels(names(seqs), seqlevelsInUse(peaks),
                                      argName="regions")]
  
  if(is.null(randomAcc))
    randomAcc <- length(regions)<1000 & maxFragLength<=10000

  if(!randomAcc){
    param <- .getBamChunkParams(bam_files[[1]], flgs=flgs,
                                keepSeqLvls=names(seqs),  nChunks=splitByChr)
  }else{
    # resize and merge regions for random access to avoid double-counting
    sizeExt <- sum(abs(c(extend+shift+maxFragLength,1L)))
    regions2 <- reduce(resize(regions, width(regions)+sizeExt, fix = "center"))
    param <- list(x=ScanBamParam(flag=flgs, which=regions2))
  }
  
  cnts <- lapply(bam_files, FUN=function(bamfile){
    if(verbose) message("Reading file ",bamfile)
    res <- lapply(names(param), FUN=function(x){
      r <- .bam2bwGetReads(bamfile, paired=paired, param=param[[x]], type=type,
                          extend=extend, shift=shift, minFragL=minFragLength,
                          maxFragL=maxFragLength, strandMode=strandMode, si=seqs)
      list(ov=countOverlaps(regions, r, type=ov.type, maxgap=maxgap,
                            ignore.strand=ignore.strand, minoverlap=minoverlap),
           reads=metadata(r)$reads)
    })
    depth <- sum(vapply(res, FUN.VALUE=integer(1), FUN=function(x) x$reads))
    res <- rowSums(vapply(res, \(x) x[[1]], integer(length(regions))))
    gc(full=TRUE, verbose=FALSE)
    list(ov=res, depth=depth)
  })
  
  depths <- vapply(cnts, FUN.VALUE=integer(1), FUN=function(x) x$depth)
  
  cnts <- matrix(unlist(lapply(cnts, \(x) x[[1]])), ncol=length(bam_files))
  if(is.null(names(bam_files))){
    if(!any(duplicated(bn<-basename(bam_files)))){
      colnames(cnts) <- gsub("\\.bam$","",bn,ignore.case=TRUE)
    }else if(!any(duplicated(bn<-dirname(bam_files)))){
      colnames(cnts) <- bn
    }else{
      colnames(cnts) <- bam_files
    }
  }else{
    colnames(cnts) <- names(bam_files)
  }
  
  row.names(cnts) <- as.character(granges(regions))
  se <- SummarizedExperiment(list(counts=cnts), rowRanges=regions)
  se$depth <- depths
  se
}


#' peakPbCountsSE
#' 
#' Generate a pseudobulk peak counts SummarizedExperiment.
#'
#' @param fragfile The path to a Tabix-indexed fragment file.
#' @param peaks A GRanges of the regions in which to count.
#' @param bcmap A named vector, indicating the pseudobulk sample (values) in 
#'   which to include each barcode (names).
#' @param genome A optional genome object or path to a genom fasta file. If
#'   included, GC bias will be added to the rowData of the output object.
#' @param insertions If TRUE, (shifted) Tn5 insertions events are counted 
#'   instead of fragments. This means that each fragment gets counted twice
#'   (for both ends). Default FALSE.
#'
#' @returns A \link[SummarizedExperiment]{RangedSummarizedExperiment} with a 
#'   'counts' assay, and columns corresponding to each unique value of `bcmap`.
#' @export
peakPbCountsSE <- function(fragfile, peaks, bcmap, genome=NULL,
                           insertions=FALSE){
  if(is.data.frame(bcmap) && !is.null(row.names(bcmap)) && 
     "group" %in% colnames(bcmap))
    bcmap <- setNames(bcmap$group, row.names(bcmap))
  
  stopifnot(length(bcmap)>1 && !is.null(names(bcmap)) &&
              (is.character(bcmap) || is.factor(bcmap)))
  
  frags <- Rsamtools::TabixFile(fragfile)
  resl <- tabixChrApply(frags, fn=function(x){
    x <- x[which(x$name %in% names(bcmap))]
    x$name <- factor(x$name, names(bcmap))
    x <- x[!is.na(x$name)]
    sapply(split(x, bcmap[as.integer(x$name)]), \(y){
      if(insertions){
        y <- epiwraps:::.align2cuts(resize(y, fix="center", width(y)-8L))
      }
      countOverlaps(peaks, y)
    })
  })
  gc(full=TRUE, verbose=FALSE)
  mat <- sapply(setNames(unique(bcmap),unique(bcmap)), \(x){
    y <- lapply(resl, \(y){
      if(!is.null(dim(y)) && x %in% colnames(y)) return(y[,x])
      NULL
    })
    y <- y[!sapply(y,is.null)]
    Reduce("+",y)
  })
  
  se <- SummarizedExperiment(list(counts=mat), rowRanges=peaks)

  if(!is.null(genome)){
    if(is.character(genome) && length(genome)==1)
      genome <- Rsamtools::FaFile(genome)
    se <- tryCatch({
      chromVAR::addGCBias(se, genome=genome)
    }, error=function(e){
      warning("Failed to add GC bias to object: ", e)
      se
    })
  }
  
  se
}