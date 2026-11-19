#' importFragments
#' 
#' A fast, lightweight replacement for \code{rtracklayer::import(format="bed")}
#' specialized for tabix-indexed fragment files (chr, start, end, barcode[, 
#' count]).
#'
#' @param x The path to a tabix-indexed fragment file, or a TabixFile object.
#' @param which A \code{\link[GenomicRanges]{GRanges}} of the region(s) for 
#'   which to read fragments. If missing, the whole file is read.
#'
#' @param useCLI Logical; whether to stream the data through the `tabix` binary
#'   (much faster, if available) rather than `rtracklayer::import`.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} with the same conventions as 
#'   \code{rtracklayer::import} for BED files (1-based, closed coordinates; 
#'   unstranded), with the barcode in the `name` column and, if present in the 
#'   file, the read count in the `score` column.
#' @export
#' @importFrom Rsamtools TabixFile
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges seqnames start end
#' @importFrom IRanges IRanges
#' @examples
#' frags <- tempfile(fileext = ".tsv")
#' d <- data.frame(chr=rep(letters[1:2], each=10), start=rep(100*(1:10),2))
#' d$end <- d$start + 15L
#' d$cell <- paste0("barcode",sample.int(3, nrow(d), replace=TRUE))
#' write.table(d, frags, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
#' frags <- Rsamtools::bgzip(frags)
#' Rsamtools::indexTabix(frags, format = "bed")
#' importFragments(frags, GRanges("a", IRanges(1, 1000)))
importFragments <- function(x, which=NULL, useCLI=nzchar(Sys.which("tabix"))){
  x <- TabixFile(x)
  d <- NULL
  if(useCLI){
    # streaming directly into fread avoids creating R strings, which is the
    # slowest step of scanTabix; this requires the `tabix` binary (htslib)
    regs <- if(is.null(which)) "" else 
      paste(shQuote(sprintf("%s:%d-%d", as.character(seqnames(which)),
                            start(which), end(which))), collapse=" ")
    d <- tryCatch(
      suppressWarnings(fread(cmd=paste("tabix", shQuote(path(x)), regs), sep="\t", 
            header=FALSE, showProgress=FALSE, fill=TRUE)),
      error=function(e) NULL)
  }
  if(is.null(d)){
    # fallback (no `tabix` binary): rtracklayer
    if(is.null(which)) return(rtracklayer::import(path(x), format="bed"))
    return(rtracklayer::import(path(x), format="bed", which=which))
  }
  if(nrow(d)==0L){
    gr <- GRanges()
    gr$name <- character()
    return(gr)
  }
  gr <- GRanges(as.character(d[[1L]]), 
                IRanges(as.integer(d[[2L]])+1L, as.integer(d[[3L]])))
  if(ncol(d)>=4L) gr$name <- as.character(d[[4L]])
  if(ncol(d)>=5L) gr$score <- as.numeric(d[[5L]])
  gr
}
