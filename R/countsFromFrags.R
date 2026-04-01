#' peakCountsFromFrags
#' 
#' Creates a SummarizedExperiment of cell-level or pseudo-bulk-level fragment 
#' (or insertion) counts from a tabix-indexed fragment file and a set of 
#' regions of interest.
#'
#' @param fragfile A path to the tabix-indexed fragment file.
#' @param regions A \code{\link[GenomicRanges]{GRanges}} of regions in which to 
#'   count overlaps.
#' @param barcodes An optional character vector of cell barcodes to include. 
#'   If provided, only these barcodes will be considered, if `NULL`, all 
#'   barcodes included in the file are used. If `barcodes` is a named vector,
#'   the names will be considered to represent the cell barcode, the values
#'   to represent the pseudo-bulk sample in which to include the respective 
#'   barcodes, and pseudo-bulk counts will be returned.
#' @param insertions Logical; if \code{TRUE}, the ends of the fragments 
#'   (insertions) are counted instead of the entire fragment. This means each 
#'   fragment can contribute up to two counts. Default \code{FALSE}.
#' @param minFragLength Minimum fragment length to be considered. Default 1.
#' @param maxFragLength Maximum fragment length to be considered. Default 5000.
#' @param ... Passed to \code{\link{tabixChrApply}}.
#' @inheritParams peakCountsFromBAM
#' 
#' @return A \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}} 
#'   with a 'counts' assay. The mean fragment length per region is also stored
#'   in the `rowData` of the object.
#' @export
#' @importFrom GenomicRanges GRanges findOverlaps width granges
#' @importFrom IRanges IRanges ranges
#' @importFrom S4Vectors queryHits subjectHits splitAsList
#' @importFrom Matrix sparseMatrix
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom stats median
#' @examples
#' # generate dummy regions and save them to a temp file:
#' frags <- tempfile(fileext = ".tsv")
#' d <- data.frame(chr=rep(letters[1:2], each=10), start=rep(100*(1:10),2))
#' d$end <- d$start + 15L
#' d$cell <- paste0("barcode",sample.int(3, nrow(d), replace=TRUE))
#' write.table(d, frags, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
#' # tabix-index it
#' frags <- Rsamtools::bgzip(frags)
#' Rsamtools::indexTabix(frags, format = "bed")
#' # we create regions of interest:
#' regions <- GRanges(c("a","b"), IRanges(400,width=300))
#' # we get the counts:
#' se <- peakCountsFromFrags(frags, regions)
#' se
#' # we could also get pseudobulk counts by passing a barcode map:
#' bcmap <- setNames(c("PB1","PB1","PB2"),paste0("barcode",1:3))
#' pb <- peakCountsFromFrags(frags, regions, barcodes=bcmap)
#' pb
peakCountsFromFrags <- function(fragfile, regions, barcodes = NULL, 
                                insertions = FALSE, minFragLength = 1L, 
                                maxFragLength = 5000L, ov.type="any",
                                maxgap=-1L, minoverlap=1L,
                                ignore.strand=TRUE, ...){
  
  stopifnot(is(regions, "GRanges"))
  if(!is(fragfile, "TabixFile")){
    if(length(fragfile)>1)
      stop("Please pass a single filepath into `fragfile`.")
    stopifnot(file.exists(fragfile))
  }

  if(!is.null(barcodes)){
    if(!is.character(barcodes)) barcodes <- as.character(barcodes)
    if(!is.null(names(barcodes))){
      stopifnot(sum(duplicated(names(barcodes)))==0)
    }else{
      barcodes <- unique(barcodes)
      names(barcodes) <- barcodes
    }
  }
  bc_lvls <- unique(barcodes)
  
  depth <- NULL
  
  chunk_fn <- function(x) {
    wi <- width(x)
    toKeep <- which(wi >= minFragLength & wi <= maxFragLength)
    if(!is.null(barcodes)){
      toKeep <- intersect(toKeep, which(x$name %in% names(barcodes)))
      x <- x[toKeep]
      x$name <- factor(x$name, names(barcodes), as.character(barcodes))
      factor(1:5, 1:10, c(rep(1:3,each=3),8))
      depth <- as.integer(table(x$name))
    }else{
      x <- x[toKeep]
    }
    wi <- wi[toKeep]
    
    if (insertions) x <- .align2cuts(x)
    
    ov <- findOverlaps(regions, x, maxgap=maxgap, minoverlap=minoverlap,
                       ignore.strand=ignore.strand)
    
    co <- data.frame(i=queryHits(ov))
    if(!is.null(barcodes)){
      co$bc <- as.integer(x$name[subjectHits(ov)])
    }else{
      co$bc <- as.character(x$name[subjectHits(ov)])
    }
    fl_vec <- rep(0L, length(regions))
    if(length(ov) > 0){
      fl <- sum(splitAsList(wi[subjectHits(ov)],queryHits(ov)))
      fl_vec[unique(sort(queryHits(ov)))] <- fl
    }
    return(list(counts=co, fl=fl_vec, depth=depth))
  }

  res <- tabixChrApply(fragfile, fn=chunk_fn, ...)

  if(!is.null(barcodes)) depth <- Reduce("+", lapply(res, \(x) x$depth))
  fl <- Reduce("+",lapply(res, \(x) x$fl))
  res <- dplyr::bind_rows(lapply(res, \(x) x$counts))
  if(is.null(res) || length(res)==0){
    if(is.null(barcodes)){
      mat <- sparseMatrix(i=integer(0), j=integer(0),
                          dims=c(length(regions), 1))
    }else{
      mat <- sparseMatrix(i=integer(0), j=integer(0), 
                          dims=c(length(regions), length(bc_lvls)),
                          dimnames=list(NULL, bc_lvls))
    }
  }else{
    if(is.null(barcodes)){
      bc_lvls <- barcodes <- unique(sort(res$bc))
      res$bc <- factor(res$bc, barcodes)
    }
    mat <- sparseMatrix(i=res$i, j=as.integer(res$bc), x=1,
                        dims=c(length(regions), length(bc_lvls)),
                        dimnames=list(NULL, bc_lvls))
    fl <- fl/Matrix::rowSums(mat)
    regions$flbias <- log10(1+fl)
  }
  row.names(mat) <- as.character(granges(regions))
  se <- SummarizedExperiment(list(counts=mat), rowRanges=regions)
  if(!is.null(depth)) se$depth <- depth
  se
}