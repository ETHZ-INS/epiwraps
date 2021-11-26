#' annotateRegions
#' 
#' Annotates a GRanges on the basis of an annotation object (e.g. 
#' \code{\link[ensembldb]{EnsDb}}).
#' 
#' @param regions A GRanges object
#' @param anno An annotation object, such as an \code{\link[ensembldb]{EnsDb}})
#' object, a TxDb object, or a GRanges object of a GENCODE-like gtf or the 
#' path to such a file.
#' @param proximal The threshold(s) for TSS proximal regions. Multiple values
#' will result in multiple class factor levels.
#' @param filter An \code{\link[AnnotationFilter]{AnnotationFilter}} to filter
#' transcripts. Only used if `anno` is an \code{\link[ensembldb]{EnsDb}}).
#' @param extra An optional named list of GRanges for additional overlaps. Each
#' list element will create an additional binary metadata column.
#' @param ignore.strand Whether to ignore the strand for the overlap with 
#' elements of `extra` (default TRUE).
#' @param ... Passed to \code{\link[GenomicRanges]{overlapsAny}} for the 
#' overlaps with `extra`.
#'
#' @return The sorted `regions` object with additional annotation columns.
#' 
#' @import GenomicRanges
#' @importFrom GenomeInfoDb keepSeqlevels seqlevelsStyle
#' @importFrom ensembldb transcripts exons 
#' @importFrom rtracklayer import import.gff
#' @importFrom AnnotationFilter AnnotationFilterList
#' @export
annotateRegions <- function(regions, anno, proximal=c(2500,1000), 
                            filter=AnnotationFilterList(),
                            extra=list(), ignore.strand=TRUE, ...){
  stopifnot(is(regions, "GRanges"))
  stopifnot(is.list(extra))
  if(length(extra)>0 && is.null(names(list())))
    stop("`extra` should be a named list.")
  extra <- lapply(extra, FUN=function(x){
    if(is.character(x)) x <- rtracklayer::import(x)
    if(is(x,"GRangesList")) x <- unlist(x)
    if(!is(x,"GRanges")) stop("Some elements of `extra` are not GRanges or ",
                              "files containing genomic intervals.")
    x
  })
  regions <- sort(regions)
  if(is.character(anno) && length(anno)==1)
    anno <- rtracklayer::import.gff(anno)
  if(is(anno,"TxDb") || is(anno,"EnsDb")){
    if(is(anno,"TxDb")){
      tx <- transcripts(anno, columns=c("tx_name","gene_id"))
      tx$transcript_id <- tx$tx_name
      tx$tx_name <- NULL
      tx$gene_name <- tx$gene_id
      ex <- exons(anno, columns=NULL)
    }else{
      tx <- transcripts(anno, columns=c("tx_id", "gene_id", "gene_name"), 
                        filter=filter)
      tx$transcript_id <- tx$tx_id
      tx$tx_id <- NULL
      ex <- exons(anno, columns=NULL, filter=filter)
    }
    tx$type <- "transcript"
    mcols(ex)$type <- "exon"
    anno <- sort(c(tx,ex))
    anno$type <- as.factor(anno$type)
  }
  if(!is(anno, "GRanges")) stop("Unknown `anno` format!")
  seqlevelsStyle(anno) <- seqlevelsStyle(regions)
  anno <- anno[anno$type %in% c("transcript","exon","dispersed_repeat")]
  sl <- intersect(seqlevels(anno),seqlevels(regions))
  if(length(sl)==0) stop("No seqlevel in common!")
  anno <- keepSeqlevels(anno, sl, pruning.mode="coarse")
  regions <- keepSeqlevels(regions, sl, pruning.mode="coarse")
  anno$type <- relevel(droplevels(anno$type),"exon")
  tss <- anno[anno$type=="transcript"]
  tss1 <- tss[strand(tss)=="+"]
  tss2 <- tss[strand(tss)=="-"]
  tss <- c( GRanges(seqnames(tss1), IRanges(start(tss1), width=1), 
                    strand=strand(tss1), 
                    mcols(tss1)[,c("transcript_id","gene_id","gene_name")]),
            GRanges(seqnames(tss2), IRanges(end(tss2), width=1),
                    strand=strand(tss2), 
                    mcols(tss2)[,c("transcript_id","gene_id","gene_name")]) )
  d <- distanceToNearest(regions, tss)
  regions$distance2nearestTSS <- NA_integer_
  regions$distance2nearestTSS[d@from] <- mcols(d)$distance
  tmp <- cbind(start(regions)[d@from], end(regions)[d@from])-start(tss)[d@to]
  regions$distance2nearestTSS[d@from] <- apply(tmp,1,FUN=function(x){
    if(length(unique(sign(x)))>1) return(0)
    x[which.min(abs(x))]
  })
  regions$nearestTSS <- regions$nearestTSS.gene_name <- NA_character_
  mcols(regions)[d@from,c("nearestTSS","nearestTSS.gene_name")] <- 
    mcols(tss)[d@to,c("transcript_id","gene_name")]
  o <- findOverlaps(regions, anno)
  ll <- sapply(split(o@to, o@from), FUN=function(x) sort(anno$type[x])[1])
  regions$TSS.overlap <- "intergenic"
  regions$TSS.overlap[as.numeric(names(ll))] <- 
    c("exonic","intronic","repeat")[as.numeric(ll)]
  regions$class <- regions$TSS.overlap
  proximal <- sort(proximal, decreasing=TRUE)
  for(i in seq_along(proximal)){
    lab <- paste0("proximal ",
                  ifelse(i!=length(proximal),paste0(">",proximal[i+1],"&"),""),
                  "<=",proximal[i],"bp")
    regions$class[abs(regions$distance2nearestTSS)<=proximal[i]] <- lab
  }
  regions$class[abs(regions$distance2nearestTSS)==0] <- "TSS"
  regions$TSS.overlap <- factor(regions$TSS.overlap)
  regions$class <- as.factor(regions$class)
  for(x in names(extra)){
    seqlevelsStyle(extra[[x]]) <- seqlevelsStyle(regions)
    mcols(regions)[[x]] <- overlapsAny( regions, extra[[x]], 
                                        ignore.strand=ignore.strand, ... )
  }
  regions
}
