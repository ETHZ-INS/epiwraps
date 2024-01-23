#' signal2Matrix
#' 
#' Reads the signals around (the centers of) a set of regions.
#'
#' @param filepaths A named vector of filepaths (e.g. to bigwig files; bam files
#'   are also supported, but with limited functionalities). Can also be a named
#'   list including a combination of paths to such files and `GRanges` object.
#'   For `GRanges` objects, the `score` column will be used (absolute coverage
#'   mode).
#' @param regions A `GRanges` of the regions/positions around which to plot, or
#' the path to a bed file of such regions.
#' @param extend Number of basepair to extend on either side of the regions. 
#'   Must be a multiple of `w`. Can also be an integer of length 2, indicating 
#'   the extension upstream and downstream.
#' @param w Bin size
#' @param scaledBins The number of bins for the scale region (ignored if 
#'   `type="center"`) 
#' @param type Either 'center' (plots fixed-size region around the centers of 
#'   `regions`) or 'scaled' (scales the signal in `regions` and plot 
#'   surroundings)
#' @param binMethod Whether to compute the 'max' (default), 'mean' or 'min' per
#'   bin.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object, or 
#' the number of threads to use to read and prepare the data. Note that the 
#'   rate-limiting process is reading from disk, so unless you have an unusually
#'   fast disk, using multi-threading is actually likely to slow down rather 
#'   than speed up the process.
#' @param verbose Logical; whether to print processing information
#' @param ... Passed to \code{\link[EnrichedHeatmap]{as.normalizedMatrix}} when
#'   reading bigwig files, or to \code{\link{getBinSignalFromBam}} when reading
#'   bam files.
#'
#' @return A list of `normalizeToMatrix` objects
#' @export
#' 
#' @import GenomicRanges
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam
#' @importFrom Rsamtools scanBamFlag ScanBamParam countBam
#' @import EnrichedHeatmap
#' @importFrom rtracklayer import BigWigSelection
#' @importFrom GenomicAlignments readGAlignmentPairs
#' 
#' @examples 
#' # we fetch the path to the example bigwig file:
#' (bw <- system.file("extdata/example_atac.bw", package="epiwraps"))
#' # we load example regions:
#' regions <- rtracklayer::import(system.file("extdata/example_peaks.bed", 
#'                                            package="epiwraps"))
#' length(regions)
#' # we obtain the matrix of the signal around the regions:
#' m <- signal2Matrix(bw, regions)
#' dim(m[[1]])
#' # we can plot it with:
#' plotEnrichedHeatmaps(m)
#' # we could also take a broader range around the center of the regions, and 
#' # use bigger bins:
#' m <- signal2Matrix(bw, regions, extend=2000, w=20)
#' # the matrix has the same size, but shows broader regions:
#' dim(m[[1]])
#' plotEnrichedHeatmaps(m)
signal2Matrix <- function(filepaths, regions, extend=2000, w=NULL,
                          scaledBins=50L, type=c("center","scaled"),
                          binMethod=c("max","mean","min"), BPPARAM=1L, 
                          ret=c("list","ESE"), verbose=TRUE, ...){
  type <- match.arg(type)
  ret <- match.arg(ret)
  binMethod <- match.arg(binMethod)
  
  if(!all(unlist(lapply(filepaths, FUN=function(x){
    is(x, "GRanges") || file.exists(x) }))))
    stop("Some of the files given do not exist, check the paths provided.")
  if(type!="center" && any(unlist(lapply(filepaths, FUN=function(x){
    !is(x, "GRanges") && grepl("\\.bam$",x,ignore.case=TRUE)}))))
    stop("Only `type='center'` can be used for BAM files.")
  
  if(is.null(w) || is.na(w)) w <- round(max(1,mean(extend)/100))
  w <- as.integer(w)
  stopifnot(w>0)
  if(!all((extend %% w)==0)) stop("`extend` should be a multiple of `w`!")
  extend <- as.integer(extend)
  if(length(extend)==1) extend <- c(extend,extend)
  
  if(is.null(names(filepaths)))
    names(filepaths) <- .cleanFileNames(filepaths)
  
  if(is.character(regions)){
    stopifnot(is.character(regions) && length(regions)==1)
    if(!file.exists(regions)) stop("The `regions` file appears not to exist.")
    regions <- rtracklayer::import(regions)
  }
  stopifnot(is(regions,"GRanges") || is(regions, "GRangesList"))
  
  # give names to regions
  if(is.null(names(regions))){
    if(is(regions, "GRanges")){
      names(regions) <- paste0(as.character(granges(regions)))
    }else{
      names(regions) <- paste0("region", seq_along(regions))
    }
  }
  rnames <- names(regions)
  regions <- sort(regions)

  if(is(regions, "GRangesList")){
    if(type=="center") stop("A GRangesList cannot be used with type='center'.")
    regions2 <- regions
    regions <- .grlBounds(regions)
  }else if(type=="center"){
    regions2 <- resize(regions, fix="center", width=sum(extend))
    regions2 <- shift(regions2, round((extend[2]-extend[1])/2))
  }else{
    regions2 <- regions
  }
  
  if(is(regions2, "GRanges")){
    # ensure that the regions are on represented seqnames
    regions <- .checkRegions(filepaths, regions, verbose=verbose,
                              trimOOR=FALSE)
    regions2 <- regions2[names(regions2)]
  }

  if(type=="scale"){
    upstream <- flank(regions, extend[[1]], start=TRUE)
    downstream <- flank(regions, extend[[2]], start=FALSE)
  }
  
  ff <- function(filename, ...){
             
    filepath <- filepaths[[filename]]
    if(is(filepath, "GRanges")){
      if(verbose) message("Computing signal from GRanges '", filename, "'...")
    }else{
      if(verbose) message("Reading ", filepath)
      if(grepl("\\.rds$",filepath,ignore.case=TRUE))
        filepath <- readRDS(filepath)
    }
    
    if(is(filepath, "GRanges")){
      
      ####### GRanges INPUT
      if(type=="scale"){
        target_ratio <- w*scaledBins/sum(extend)
        mat <- normalizeToMatrix(filepath, regions, w=w, extend=extend, 
                                 value_column="score", mean_mode="absolute", 
                                 target_ratio=target_ratio, ...)
      }else{
        mat <- normalizeToMatrix(filepath, resize(regions,1L,fix="center"), w=w,
                                 extend=extend, value_column="score", ..., 
                                 mean_mode="absolute")        
      }
      ####### END GRanges INPUT
      
    }else if(grepl("\\.bam$",filepath,ignore.case=TRUE)){
      
      ####### BAM INPUT
      
      mat <- getBinSignalFromBam(filepath, regions2, ...)
      if(w!=1){
        nwin <- round(sum(extend)/w)
        if(nwin %% 2 != 0) nwin <- nwin+1
        mat <- sapply(split(seq_len(ncol(mat)),cut(seq_len(ncol(mat)),nwin)), 
                      FUN=function(x) rowMeans(mat[,x]) )
        extend <- floor(ncol(mat)/2)
      }
      mat <- EnrichedHeatmap::as.normalizedMatrix( 
        unclass(mat), extend=extend, signal_name=filename, k_target=0, 
        k_upstream=extend[1]/w, k_downstream=extend[2]/w+(w==1))
      
      ####### END BAM INPUT
      
    }else if(grepl("\\.bw$|\\.bigwig$", filepath, ignore.case=TRUE)){

      ####### BIGWIG INPUT
      
      if(type=="scaled"){
        co <- rtracklayer::import(filepath, format="BigWig", 
                                  selection=BigWigSelection(reduce(regions)))
        co <- coverage(co, weight=co$score)
        upstream <- .getBinSignalFromBW(co, upstream, w=w, 
                                        method=binMethod, verbose=verbose)
        downstream <- .getBinSignalFromBW(co, downstream, w=w, 
                                          method=binMethod, verbose=verbose)
        mat <- .getScaledSignalFromBW(co, regions, nBins=scaledBins, 
                                      method=binMethod, verbose=verbose)
        mat <- cbind(upstream, mat, downstream)
        mat <- EnrichedHeatmap::as.normalizedMatrix(unclass(mat), extend=extend,
                  signal_name=filename, k_target=scaledBins, ...,
                  k_upstream=extend[[1]]/w, k_downstream=extend[[2]]/w )
      }else{
        mat <- .getBinSignalFromBW(filepath, regions2, w=w, 
                                   method=binMethod, verbose=verbose)
        mat <- EnrichedHeatmap::as.normalizedMatrix(unclass(mat), extend=extend,
                  signal_name=filename, k_target=0L, ...,
                  k_upstream=extend[[1]]/w, k_downstream=extend[[2]]/w )
      }

      ####### END BIGWIG INPUT
      
    }else{
      stop("Unknown file format")
    }
    mat <- tryCatch({
        if(is.null(names(mat))) names(mat) <- names(regions)
        mat[rnames,]
      }, error=function(e){
        if(verbose) warning(e)
        mat
      })
    mat
  }

  BPPARAM <- .getBP(BPPARAM)
  if(BiocParallel::bpnworkers(BPPARAM)==1){
    ml <- lapply(setNames(names(filepaths),names(filepaths)), FUN=ff, ...)
  }else{
    ml <- bplapply(setNames(names(filepaths),names(filepaths)), FUN=ff, 
                   BPPARAM=BPPARAM, ...)
  }

  ml <- tryCatch(.comparableMatrices(ml), error=function(e){
    if(verbose) warning(e)
    ml
  })
  if(ret=="ESE"){
    ml <- tryCatch(.ml2ese(ml, rowRanges=regions), error=function(e){
      if(verbose)
        warning("Could not create ESE object (list returned instead):", e)
      ml
    })
  }
  ml
}

.filterRegions <- function(regions, seqlvls, verbose=TRUE){
  if(length(missing <- setdiff(seqlevelsInUse(regions), seqlvls))>0){
    toRemove <- seqnames(regions) %in% missing
    if(verbose) warning(length(missing),
                        " seqlevel(s) missing from the bigwig file.\n",
                        sum(toRemove), " regions on these sequences ",
                        "will be ignored.")
    regions <- regions[!toRemove]
  }
  keepSeqlevels(regions, seqlevelsInUse(regions), 
                pruning.mode="coarse")
}

.getBinSignalFromBW <- function(filepath, regions2, w=1L, 
                                method=c("max","min","mean"), verbose=TRUE){
  method <- match.arg(method)
  stopifnot(length(unique(width(regions2)))==1)
  if(is(filepath, "RleList")){
    sqlvls <- names(filepath)
    co <- filepath
  }else{
    sqlvls <- seqlevels(BigWigFile(filepath))
    co <- rtracklayer::import(filepath, format="BigWig", as="RleList",
                              selection=BigWigSelection(reduce(regions2)))
  }
  regions2 <- .filterRegions(regions2, sqlvls, verbose=verbose)
  if(w==1L){
    mat <- views2Matrix(Views(co[seqlevels(regions2)], regions2), padVal=0L)
  }else{
    mat <- .getBinSignal(co[seqlevels(regions2)], regions2, w=w, method=method)
  }
  wRev <- which(as.factor(strand(regions2))=="-")
  if(length(wRev)>0) 
    mat[wRev,] <- mat[wRev,seq(from=ncol(mat), to=1L),drop=FALSE]
  mat
}

.getScaledSignalFromBW <- function(filepath, regions2, nBins=50L, verbose=TRUE,
                                   method=c("max","min","mean"), useRle=FALSE){
  method <- match.arg(method)
  if(is(filepath, "RleList")){
    sqlvls <- names(filepath)
    co <- filepath
  }else{
    sqlvls <- seqlevels(BigWigFile(filepath))
    if(is(regions2, "GRangesList")){
      co <- rtracklayer::import(filepath, format="BigWig", as="RleList", 
                        selection=BigWigSelection(reduce(unlist(regions2))))
    }else{
      co <- rtracklayer::import(filepath, format="BigWig", as="RleList",
                                selection=BigWigSelection(reduce(regions2)))
    }
  }
  regions2 <- .filterRegions(regions2, sqlvls, verbose=verbose)
  if(is(regions2, "GRangesList")){
    strands <- unique(strand(regions2))
    # join exons coverages
    v <- lapply(split(granges(regions2), runValue(seqnames(regions)), drop=TRUE),
           FUN=function(r){
      v <- Views(co[[as.character(seqnames(r[[1]])[1])]], ranges(unlist(r)))
      .mergeRleListItems(RleList(v), rep(factor(names(r)), lengths(r)))
    })
    v <- unlist(v, recursive=FALSE, use.names=FALSE)
    mat <- do.call(rbind, lapply(v, FUN=function(x){
      t(sapply(x, FUN=function(x){
        x <- as.numeric(x)
        x <- rep(x,each=max(1L,ceiling(nBins/length(x))))
        x <- splitAsList(x, cut(seq_along(x), breaks=nBins, labels=FALSE))
        switch(method, max=max(x), min=min(x), mean=mean(x))
      }))
    }))
  }else{
    strands <- strand(regions2)
    mat <- .getBinSignal(co[seqlevels(regions2)], regions2, k=nBins,
                         method=method)
  }
  wRev <- which(as.factor(strands)=="-")
  if(length(wRev)>0) mat[wRev,] <- mat[wRev,seq(from=ncol(mat), to=1L)]
  mat
}

#' getBinSignalFromBam
#'
#' This is a wrapper around \code{\link[genomation]{ScoreMatrix}} to enable BAM
#' support in \code{\link{signal2Matrix}}.
#'
#' @param filepath The path to the (indexed) bam file
#' @param regions A `GRanges` of the regions/positions around which to plot
#' @param cuts Whether to count cuts (e.g. beginning/end of fragments) rather
#' than coverage (ignored unless the input are bam files)
#' @param RPM Whether to perform RPM normalization (for bam input)
#' @param paired Whether to consider whole fragments
#' @param flgs Flags for bam reading
#'
#' @return A matrix
#' @export
getBinSignalFromBam <- function(filepath, regions, cuts=FALSE, RPM=TRUE, 
                                paired=TRUE, ...,
                                flgs=scanBamFlag(isDuplicate=FALSE,
                                                 isSecondaryAlignment=FALSE)){
  if(!suppressWarnings(requireNamespace("genomation", quietly=TRUE)))
    stop("Please install the 'genomation' package to enable plotting 
             directly from bam files.")
  libsize <- NULL
  if(is.null(RPM) || RPM)
    libsize <- Rsamtools::countBam(filepath, 
                                   param=ScanBamParam(flag=flgs))$records
  if(cuts){
    params <- Rsamtools::ScanBamParam(which=regions, flag=flgs)
    bam <- GenomicAlignments::readGAlignmentPairs(filepath, param=params)
    bam <- coverage(.align2cuts(bam))
    
    mat <- genomation:::ScoreMatrix(bam, windows=regions, ...,
                                   strand.aware=FALSE, library.size=libsize)
    rm(bam)
    if(is.null(RPM) || RPM) mat <- mat*1000000/libsize
  }else{
    mat <- genomation:::ScoreMatrix(filepath, regions, bam.paired.end=paired, 
                                   unique=TRUE, library.size=libsize, ...)
    if(is.null(RPM) || RPM) mat <- mat*1000000/libsize
  }
  mat
}

#' @importFrom EnrichedHeatmap makeWindows
.getBinSignal <- function(co, regions, w=NULL, k=NULL, short.keep=TRUE,
                          method=c("max","mean","min")){
  method <- match.arg(method) 
  if(is.null(names(regions))) names(regions) <- as.character(granges(regions))
  norder <- names(regions)
  regions <- sort(regions)
  v <- Views(co, makeWindows(regions, w=w, k=k, short.keep=short.keep))
  v <- switch(method, min=viewMins(v), max=viewMaxs(v), mean=viewMeans(v))
  v <- v[lengths(v)>0]
  v <- unlist(v)
  if(all(w <- is.infinite(v))){
    finiteMin <- 0L
  }else{
    finiteMin <- min(min(v[!w]),0L)
  }
  v[w] <- finiteMin
  mat <- matrix(unlist(v), nrow=length(regions), byrow=TRUE,
                dimnames=list(names(regions),NULL))
  mat[norder,]
}

# boundaries of each element of a GRangesList (e.g. transcript coords)
.grlBounds <- function(regions){
  if(!all(unique(lengths(strands <- unique(strand(regions))))==1L) |
     !all(unique(lengths(seqlvls <- unique(seqnames(regions))))==1L))
    stop("Some elements of 'regions' contain GRanges from more than one ",
         "strand and/or seqlevel.")
  x <- GRanges(seqlvls, strand=strands,
               ranges=IRanges(min(start(regions)), max(end(regions))))
  mcols(x) <- mcols(regions)
  x
}

# concatenates elements of a rlelist
.mergeRleListItems <- function(rlelist, by){
  rvs <- runValue(rlelist)
  rls <- runLength(rlelist)
  rvs <- splitAsList(unlist(rvs), rep(by,lengths(rvs)))
  rls <- splitAsList(unlist(rls), rep(by,lengths(rls)))
  RleList(lapply(setNames(seq_along(rvs),names(rvs)), FUN=function(i){
    Rle(rvs[[i]], rls[[i]])
  }))
}
