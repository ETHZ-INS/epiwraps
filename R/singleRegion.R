#' plotSignalTracks
#' 
#' A wrapper around `Gviz` for quick plotting of genomic signals in a single
#' region.
#'
#' @param files A named list or vector of paths to signal files  (e.g. 
#' bigwig/bam, but also bed files). If a list, list elements will be overlaid 
#' or aggregated (depending on the `aggregation` argument). Objects accepted by
#' \code{\link[Gviz]{DataTrack}}'s `range` argument are also accepted.
#' @param region A genomic region, either as a `GRanges` object or as a string 
#' (i.e. `region="chr5:10000-12000`). Alternatively, if `ensdb` is provided, a 
#' gene name can be given, and the gene's coordinates will be used as region.
#' @param colors Signal color(s); will be recycled for elements of `files`
#' @param type Signal plot type(s); will be recycled for elements of `files`
#' @param overlay.alpha Transparency (0 to 250) when overlaying tracks.
#' @param ensdb An optional \code{\link[ensembldb]{EnsDb}} object form which 
#' to grab transcripts.
#' @param genomeAxis Whether to plot a genome axis. Alternatively, a numeric 
#' scalar between 0 and 1 can be given, in which case a scale will be 
#' plotted of this relative size.
#' @param extend Either an integer or vector of two integers indicating the 
#' number of base pairs by which to extent on either side.
#' @param aggregation Method for aggregation data tracks, one of: 'mean' 
#' (default), 'median', 'max', 'overlay' or 'heatmap'.
#' @param transcripts Whether to show transcripts (reguires `ensdb`) as "full",
#' "collapsed" (default), "coding" (only coding transcripts) or "none".
#' Alternatively, can be a custom \code{\link[Gviz](GeneRegionTrack)} object.
#' @param genes.params Named list of parameters passed to 
#' \code{\link[Gviz](GeneRegionTrack)}.
#' @param tracks.params Named list of parameters passed to 
#' \code{\link[Gviz](DataTrack)}.
#' @param extraTracks List of extra custom tracks to be plotted.
#' @param background.title The background color of the track titles.
#' @param col.axis The color of the axes.
#' @param col.title The color of the track titles.
#' @param cex.title Expension factor for the font size of the track titles.
#' @param bed.rotation.title Rotation for track titles of bed files.
#' @param ... Passed to \code{\link[Gviz](plotTracks)}.
#'
#' @return A list of GenomeGraph tracks to be plotted.
#' 
#' @importFrom GenomicRanges reduce seqnames start end
#' @importFrom ensembldb getGeneRegionTrackForGviz genes
#' @importFrom S4Vectors mcols
#' @importFrom ensembldb getGeneRegionTrackForGviz
#' @importFrom Gviz plotTracks DataTrack OverlayTrack GeneRegionTrack 
#' @importFrom Gviz GenomeAxisTrack AnnotationTrack
#' @importFrom matrixStats rowMins rowMaxs rowMedians
#'
#' @export
plotSignalTracks <- function(files, region, colors="darkblue", ensdb=NULL, 
                             type="h",  genomeAxis=0.3, extend=0,
                             aggregation=c("mean","median", "sum", "max", 
                                           "min", "heatmap", "overlay"),
                             transcripts=c("collapsed","full","coding","none"), 
                             genes.params=list(col.line="grey40", col=NULL, 
                                               fill="#000000"),
                             tracks.params=list(), extraTracks=list(), 
                             background.title="white", col.axis="grey40", 
                             bed.rotation.title=0, col.title="black", 
                             cex.title=0.65, overlay.alpha=100, ...){
  options(ucscChromosomeNames=FALSE)
  if(!is.function(aggregation)) aggregation <- match.arg(aggregation)
  if(!is(transcripts, "GeneRegionTrack"))
    transcripts <- match.arg(transcripts)
  if(!is.null(ensdb)) stopifnot(is(ensdb,"EnsDb"))
  
  # region of interest
  region <- .parseRegion(region, ensdb)
  
  # check file formats (from names)
  fm <- lapply(files, .parseFiletypeFromName)
  if(any(lengths(lapply(fm,unique))!=1))
    stop("Cannot aggregate files of different formats!")
  
  # creating names if not specified
  if(is.null(names(files))){
    if(is.character(files)){
      names(files) <- gsub("\\.bigwig$|\\.bw$|\\.bam$|\\.bed$", "", 
                           basename(files), ignore.case=TRUE)
    }else if(is.list(files)){
      stopifnot(all(unlist(lapply(files, is.character))))
      names(files) <- sapply(seq_along(files), FUN=function(x){
        if(length(files[[x]])==1)
          return(gsub("\\.bigwig$|\\.bw$|\\.bam$|\\.bed$", "", 
                      files[[x]], ignore.case=TRUE))
        paste0("track",x)
      })
    }else{
      stop("Invalid `files` argument!")
    }
  }
  if(!is.list(files)) files <- as.list(files)
  if(any(lengths(fm)>1 && sapply(fm, FUN=function(x) any(x=="bam"))))
    warning("It is not advised to overlay/aggregate signals from bam files, ",
            "as these are not normalized.")
  
  # Handling the gene track
  gt <- NULL
  if(!is.null(ensdb) && !is(transcripts, "GeneRegionTrack") && 
     transcripts!="none"){
    ggr <- getGeneRegionTrackForGviz(ensdb, chromosome=region[[1]],
                                     start=region[[2]], end=region[[3]],
                                     featureIs=ifelse(
                                       transcripts=="collapsed",
                                       "gene_biotype","tx_biotype" ) )
    if(is.null(genes.params$transcriptAnnotation))
      genes.params$transcriptAnnotation <- ifelse(transcripts=="collapsed",
                                                  "symbol","transcript")
    if(transcripts=="collapsed") genes.params$collapseTranscripts="meta"
    if(transcripts=="coding")
      ggr <- ggr[ggr$transcript %in% 
                   unique(ggr$transcript[ggr$feature=="protein_coding"])]
    if(is.null(genes.params$name)) genes.params$name <- genome(ensdb)[[1]]
    genes.params$range <- ggr
    gt <- do.call(GeneRegionTrack, genes.params)
  }else if(is(transcripts, "GeneRegionTrack")){
    gt <- transcripts
  }
  
  # recycle colors and types
  if((rat <- length(files)/length(colors))>1)
    colors <- rep(colors,ceiling(rat))[seq_along(files)]
  if(!is.null(names(colors)) || !all(names(files) %in% names(colors)))
    names(colors) <- names(files)
  if((rat <- length(files)/length(type))>1)
    type <- rep(type,ceiling(rat))[seq_along(files)]
  if(!is.null(names(type)) || !all(names(files) %in% names(type)))
    names(type) <- names(files)
  
  tracks <- lapply(setNames(names(files),names(files)), FUN=function(subf){
    isMult <- length(files[[subf]])>1
    
    if(any(isBed <- grepl("\\.bed$",files[[subf]],ignore.case=TRUE))){
      if(isMult){
        stop("Bed merging not yet implemented, provide them as separate items.")
      }
      return(AnnotationTrack(files[[subf]], fill=colors[subf], col=NULL, 
                             rotation.title=bed.rotation.title,
                             name=ifelse(is.na(bed.rotation.title),"",subf)))
    }
    thecol <- ifelse(!isMult || aggregation!="overlay", colors[subf],
                     .maketrans(colors[subf],overlay.alpha))
    tp <- tracks.params
    tp$type <- type[[subf]]
    tp$stream <- TRUE
    tp$col <- thecol
    tp$name <- subf
    if(!isMult || aggregation=="overlay"){
      tr <- lapply(files[[subf]], FUN=function(x){
        tp$range <- x
        do.call(DataTrack,tp)
      })
      if(length(tr)>1)
        return(OverlayTrack(tr, name=subf, type=type[[subf]], title=subf,
                            fill=thecol, col=thecol))
      return(tr[[1]])
    }
    gr <- .signalsAcrossSamples(files[[subf]], region)
    if(aggregation=="heatmap"){
      tp$type <- "heatmap"
      tp$range <- gr
    }else{
      if(!is.function(aggregation))
        aggregation <- switch(aggregation, min=matrixStats::rowMins, 
                              max=matrixStats::rowMaxs, mean=rowMeans,
                              sum=rowSums, median=matrixStats::rowMedians)
      sig <- aggregation(as.matrix(mcols(gr)))
      mcols(gr) <- NULL
      gr$score <- sig
      tp$range <- gr
    }
    do.call(DataTrack,tp)
  })
  
  if(is.null(names(extraTracks)) && length(extraTracks)>0)
    names(extraTracks) <- paste0("extra",seq_along(extraTracks))
  extraTracks <- lapply(names(extraTracks), FUN=function(tname){
    x <- extraTracks[[tname]]
    if(is.character(x))
      x <- AnnotationTrack(x, name=ifelse(grepl("^extra[0-9]",tname),"",tname),
                           col=NULL, fill="grey30", 
                           rotation.title=bed.rotation.title)
    x
  })
  
  ga <- NULL
  if(!isFALSE(genomeAxis)){
    if(is.numeric(genomeAxis) & genomeAxis>0 & genomeAxis<=1){
      ga <- Gviz::GenomeAxisTrack(scale=genomeAxis, labelPos="below")
    }else{
      ga <- Gviz::GenomeAxisTrack()
    }
  }
  
  if(length(extend)<2) extend <- c(extend,extend)
  plotTracks(c(tracks, extraTracks, gt, ga), chromosome=region[[1]], 
             from=region[[2]]-extend[1], to=region[[3]]+extend[2], 
             background.title=background.title, col.axis=col.axis, 
             col.title=col.title, cex.title=cex.title, ...)
}



#' @importFrom IRanges IRanges
#' @import GenomicRanges
#' @importFrom S4Vectors mcols
#' @importFrom rtracklayer import.bw
.signalsAcrossSamples <- function(files, region){
  if(is.list(region))
    region <- GRanges(region[[1]], IRanges(region[[2]], region[[3]]))
  stopifnot(length(region)==1)
  gp <- GPos(seqnames(region), start(region):end(region))
  files <- lapply(files, which=region, rtracklayer::import.bw)
  grs <- lapply(files,FUN=function(x) x[x$score>0])
  gr <- disjoin(unlist(GRangesList(grs)), ignore.strand=TRUE)
  m <- sapply(grs, FUN=function(x){
    o <- findOverlaps(gr,x)
    y <- rep(0,length(gr))
    y[o@from] <- x$score[o@to]
    y
  })
  mcols(gr) <- m
  gr
}


#' @importFrom AnnotationFilter SymbolFilter GeneIdFilter TxIdFilter
.parseRegion <- function(region, ensdb=NULL){
  stopifnot(length(region)==1)
  if(is(region,"GRanges")) region <- as.character(region)
  stopifnot(is.character(region))
  region <- strsplit(gsub("-",":",region),":")[[1]]
  if(length(region)==1){
    if(is.null(ensdb))
      stop("`ensdb` is required when defining the region with a gene name.")
    region <- reduce(genes(ensdb, filter=SymbolFilter(region)))
    if(length(region)==0) 
      region <- reduce(genes(ensdb, filter=GeneIdFilter(region)))
    if(length(region)==0) 
      region <- reduce(genes(ensdb, filter=TxIdFilter(region)))
    if(length(region)==0){
      region <- reduce(genes(ensdb, filter=SymbolFilter(region,"startsWith")))
      if(length(region)>1)
        stop("Gene not found with this exact name, and mutliple genes match ",
             "this string.")
    }
    if(length(region)==0) stop("Gene not found!")
    if(length(region)>1)
      stop("Region input is ambiguous (multiple non-overlapping regions)")
    region <- strsplit(gsub("-",":",as.character(region)),":")[[1]][1:3]
  }
  stopifnot(length(region)==3)
  region <- as.list(region)
  region[[2]] <- as.integer(region[[2]])
  region[[3]] <- as.integer(region[[3]])
  region
}