#' plotSignalTracks
#' 
#' A wrapper around `Gviz` for quick plotting of genomic signals in a single
#' region.
#'
#' @param files A named list or vector of paths to signal files  (e.g. 
#' bigwig/bam, but also bed files). If a list, list elements will be overlaid 
#' or aggregated (depending on the `aggregation` argument). Formats accepted by
#' \code{\link[Gviz]{DataTrack}}'s `range` argument are also accepted. Can also
#' include `GRanges ` objects (which will be plotted as 
#' \code{\link[Gviz]{AnnotationTrack}}) or objects inheriting the 
#' \code{\link[Gviz]{GdObject}} class (i.e. any `Gviz` track object).
#' @param region A genomic region, either as a `GRanges` object or as a string 
#' (i.e. `region="chr5:10000-12000`). Alternatively, if `ensdb` is provided, a 
#' gene name can be given, and the gene's coordinates will be used as region.
#' @param colors Signal color(s); will be recycled for elements of `files`
#' @param type Signal plot type(s); will be recycled for elements of `files`.
#' This is ignored for bed-like files, which are shown as 
#' \code{\link[Gviz]{AnnotationTrack}}. See the `type` options of 
#' \code{\link[Gviz]{DataTrack}}. In addition to these options, the type 
#' 'alignments' can be given for bam files, which will display them as 
#' \code{\link[Gviz]{AlignmentsTrack}}.
#' @param overlay.alpha Transparency (0 to 250) when overlaying tracks.
#' @param ensdb An optional \code{\link[ensembldb]{EnsDb}} object form which 
#' to grab transcripts.
#' @param genomeAxis Whether to plot a genome axis. Alternatively, a numeric 
#' scalar between 0 and 1 can be given, in which case a scale will be 
#' plotted of this relative size.
#' @param extend Either an integer or vector of two integers indicating the 
#' number of base pairs by which to extent on either side. If `extend`<=1, this
#' will be interpreted as a fraction of the plotted region.
#' @param aggregation Method for aggregation data tracks, one of: 'mean' 
#' (default), 'median', 'max', 'overlay', 'heatmap', or 'heatmap+mean'. The 
#' latter will create a mean plot of type `type` followed by a heatmap.
#' @param transcripts Whether to show transcripts (reguires `ensdb`) as "full",
#' "collapsed" (default), "coding" (only coding transcripts) or "none".
#' Alternatively, can be a custom \code{\link[Gviz]{GeneRegionTrack}} object.
#' @param genes.params Named list of parameters passed to 
#' \code{\link[Gviz]{GeneRegionTrack}}.
#' @param tracks.params Named list of parameters passed to 
#' \code{\link[Gviz]{DataTrack}}.
#' @param align.params Named list of parameters passed to 
#' \code{\link[Gviz]{AlignmentsTrack}}. Only used for plotting bam files with
#' `type="alignments"`.
#' @param extraTracks List of extra custom tracks to be plotted.
#' @param background.title The background color of the track titles.
#' @param col.axis The color of the axes.
#' @param col.title The color of the track titles.
#' @param cex.title Expension factor for the font size of the track titles.
#' @param bed.rotation.title Rotation for track titles of bed files.
#' @param ... Passed to \code{\link[Gviz]{plotTracks}}.
#'
#' @return A list of GenomeGraph tracks to be plotted.
#' 
#' @importFrom GenomicRanges reduce seqnames start end
#' @importFrom ensembldb getGeneRegionTrackForGviz genes
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb genome seqlevels
#' @importFrom Gviz plotTracks DataTrack OverlayTrack GeneRegionTrack 
#' @importFrom Gviz GenomeAxisTrack AnnotationTrack AlignmentsTrack
#' @importFrom matrixStats rowMins rowMaxs rowMedians
#'
#' @export
#' @examples 
#' # fetch path to example bigwig file:
#' (bw <- system.file("extdata/example_rna.bw", package="epiwraps"))
#' plotSignalTracks(list(track1=bw), region="8:22165140-22212326")
#' # if we had an EnsDb object loaded, we could just input a gene instead of 
#' # coordinates, and the transcript models would automatically show (not run):
#' # plotSignalTracks(list(track1=bw), region="BMP1", ensdb=ensdb)
#' # show all transcript variants:
#' # plotSignalTracks(list(tracks=bw), region="BMP1", ensdb=ensdb,
#' #                  transcripts="full")
plotSignalTracks <- function(files=list(), region, ensdb=NULL, colors="darkblue",
                             type="histogram",  genomeAxis=0.3, extend=0.15,
                             aggregation=c("mean","median", "sum", "max", 
                                           "min", "heatmap", "overlay",
                                           "heatmap+mean"),
                             transcripts=c("collapsed","full","coding","none"), 
                             genes.params=list(col.line="grey40", col=NULL,
                                               fill="#000000"),
                             align.params=list(color=NULL),
                             tracks.params=list(), extraTracks=list(), 
                             background.title="white", col.axis="grey40", 
                             bed.rotation.title=0, col.title="black", 
                             cex.title=0.65, overlay.alpha=100, 
                             normFactors=NULL, ...){
  if(length(files)==0 && is.null(ensdb))
    stop("No track to plot!")
  options(ucscChromosomeNames=FALSE)
  if(!is.function(aggregation)) aggregation <- match.arg(aggregation)
  if(!is(transcripts, "GeneRegionTrack"))
    transcripts <- match.arg(transcripts)
  if(!is.null(ensdb)) stopifnot(is(ensdb,"EnsDb") || is(ensdb,"TxDb"))
  
  # region of interest
  region <- .parseRegion(region, ensdb)
  if(all(extend<=1) & all(extend>=0))
    extend <- round(extend*(region[[3]]-region[[2]]))
  if(length(extend)<2) extend <- c(extend,extend)
  region2 <- list(region[[1]], max(1L, region[[2]]-extend[1]), 
                  region[[3]]+extend[2])
  
  # check file formats (from names)
  fm <- lapply(files, .parseFiletypeFromName, grOk=TRUE, trackOk=TRUE, 
               covOk=TRUE)
  if(length(files)>0){
    if(any(lengths(lapply(fm,unique))!=1))
      stop("Cannot aggregate files of different formats!")
    
    # creating names if not specified
    if(is.null(names(files))){
      if(is.character(files)){
        names(files) <- .cleanFileNames(files)
      }else if(is.list(files)){
        stopifnot(all(unlist(lapply(files, is.character))))
        names(files) <- sapply(seq_along(files), FUN=function(x){
          if(length(files[[x]])==1)
            return(.cleanFileNames(files[[x]]))
          paste0("track",x)
        })
      }else{
        stop("Invalid `files` argument!")
      }
    }
    if(!is.list(files)) files <- as.list(files)
    if(any(lengths(fm)>1 & sapply(fm, FUN=function(x) any(x=="bam"))))
      warning("It is not advised to overlay/aggregate signals from bam files, ",
              "as these are not normalized.")
    if(any(lengths(fm)>1) && aggregation=="overlay" && !is.null(normFactors))
      warning("Custom normalization factors not yet implemented for 'overlay' aggregation.")
  }
  
  # converting RleLists to temporary bigwigs
  if(length(w <- which(fm=="cov"))>0 && 
     any(sapply(files[w], FUN=object.size)>10^6))
    message("Writing coverage objects to temporary bigwigs ",
            "(this might be suboptimal for large objects)...")
  for(i in w){
    stopifnot(!is.list(files[[i]]))
    fn <- tempfile("cov",fileext=".bw")
    rtracklayer::export.bw(files[[i]], fn)
    fm[i] <- "bw"
    files[i] <- fn
  }
  
  # Handling the gene track
  gt <- NULL
  if(!is.null(ensdb) && !is(transcripts, "GeneRegionTrack") && 
     transcripts!="none"){
    ggr <- getGeneRegionTrackForGviz(ensdb, chromosome=region[[1]],
                                     start=region2[[2]], end=region2[[3]],
                                     featureIs=ifelse(transcripts=="collapsed",
                                                      "gene_biotype","tx_biotype") )
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
    ft <- .parseFiletypeFromName(files[[subf]], grOk=TRUE, trackOk=TRUE)
    if(length(ft)>1 && any(ft=="track")) stop("Custom tracks cannot be merged.")
    if(all(ft=="track")) return(files[[subf]])
    isMult <- !is(files[[subf]], "GRanges") && length(files[[subf]])>1
    if(any(ft %in% c("bed","GRanges"))){
      if(isMult && any(ft=="bed")){
        files[[subf]] <- 
          unlist(GRangesList(lapply(files[[subf]], rtracklayer::import)))
      }
      return(AnnotationTrack(files[[subf]], fill=colors[subf], col=NULL, 
                             rotation.title=bed.rotation.title,
                             name=ifelse(is.na(bed.rotation.title),"",subf)))
    }
    thecol <- ifelse(!isMult || aggregation!="overlay", colors[subf],
                     .maketrans(colors[subf],overlay.alpha))
    tp <- tracks.params
    tp$type <- type[[subf]]
    if(is.null(tp$lwd) && tp$type=="histogram") tp$lwd <- 0
    tp$stream <- TRUE
    if(is.null(tp$col)) tp$col <- thecol
    if(is.null(tp$fill)) tp$fill <- thecol
    tp$name <- subf
    if(!isMult || aggregation=="overlay"){
      tr <- lapply(files[[subf]], FUN=function(x){
        tp$range <- x
        if(grepl("^alignment", tp$type, ignore.case=TRUE)){
          if(.parseFiletypeFromName(x)=="bam"){
            ap <- align.params
            ap$range <- x
            ap$name <- subf
            ap$fill <- colors[subf]
            return(do.call(AlignmentsTrack, ap))
          }
          warning("Alignment type can only be given for bam files.")
          tp$type <- "histogram"
        }
        do.call(DataTrack, tp)
      })
      if(length(tr)>1)
        return(OverlayTrack(tr, name=subf, type=type[[subf]], title=subf,
                            fill=thecol, col=thecol))
      return(tr[[1]])
    }
    gr <- signalsAcrossSamples(files[[subf]], region2)
    if(!is.null(normFactors)){
      if(!all(names(files[[subf]]) %in% names(normFactors))){
        warning("Tracks names not found in `normFactors` - the normalization ",
                "factors will be ignored.")
      }else{
        for(f in colnames(mcols(gr))){
          mcols(gr)[[f]] <- mcols(gr)[[f]]*normFactors[[f]]
        }
      }
    }
    if(aggregation=="heatmap+mean"){
      gr2 <- gr
      mcols(gr2) <- NULL
      score(gr2) <- rowMeans(as.matrix(mcols(gr)))
      tp2 <- tp
      tp2$range <- gr2
      tp$type <- "heatmap"
      tp$range <- gr
      out <- list()
      out[[paste0(subf, "\nmean")]] <- do.call(DataTrack, tp2)
      out[[subf]] <- do.call(DataTrack, tp)
      return(out)
    }
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
      score(gr) <- sig
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
  tmns <- as.character(unlist(lapply(tracks,names)))
  tracks <- unlist(tracks, recursive = FALSE, use.names = FALSE)
  names(tracks) <- tmns
  
  plotTracks(c(tracks, extraTracks, gt, ga), 
             chromosome=region2[[1]], from=region2[[2]], to=region2[[3]], 
             background.title=background.title, col.axis=col.axis, 
             col.title=col.title, cex.title=cex.title, ...)
}


#' signalsAcrossSamples
#' 
#' Obtain a matrix of score/coverages across a region for a list of BigWig files.
#' 
#' @param files A named list of paths to biwgig files or of `GRanges` objects
#' with a `score` column.
#' @param region The region of interest, either given as a string (in the 
#' "chr:start-end" format) or as a `GRanges` of length 1.
#' @param ignore.strand Logical; whether to merge scores from the two strands
#' given stranded objects.
#' 
#' @return A disjoined `GRanges object` with the scores as metadata columns.
#'
#' @importFrom IRanges IRanges
#' @import GenomicRanges
#' @importFrom S4Vectors mcols
#' @importFrom rtracklayer import.bw
signalsAcrossSamples <- function(files, region, ignore.strand=TRUE){
  region <- .parseRegion(region, asGR=TRUE)
  files <- lapply(files, which=region, rtracklayer::import.bw)
  grs <- lapply(files,FUN=function(x) x[x$score>0])
  grs <- lapply(grs, FUN=function(x){
    keepSeqlevels(x, seqnames(region), pruning.mode="coarse")
  })
  gr <- disjoin(unlist(GRangesList(grs)), ignore.strand=ignore.strand)
  m <- sapply(grs, FUN=function(x){
    o <- findOverlaps(gr,x, ignore.strand=ignore.strand)
    y <- rep(0,length(gr))
    y[o@from] <- x$score[o@to]
    y
  })
  mcols(gr) <- m
  gr
}


#' @importFrom AnnotationFilter SymbolFilter GeneIdFilter TxIdFilter
.parseRegion <- function(region, ensdb=NULL, asGR=FALSE){
  if(is.list(region) && length(region)==3 && all(lengths(region)==1) &&
     all(is.numeric(unlist(region[2:3])))){
    if(asGR) return(GRanges(region[[1]], IRanges(region[[2]], region[[3]])))
    return(region)
  }
  stopifnot(length(region)==1)
  if(is(region,"GRanges")){
    if(asGR) return(region)
    region <- as.character(GRanges(region))
  }
  stopifnot(is.character(region))
  region <- strsplit(gsub("-",":",region),":")[[1]]
  if(length(region)==1){
    # assumes an ID is given; check in that order:
    # gene symbols, gene ids, transcript ids, or partial gene symbol matches
    if(is.null(ensdb))
      stop("`ensdb` is required when defining the region with a gene name.")
    gname <- region
    if(is(ensdb, "EnsDb")){
      region <- reduce(genes(ensdb, filter=SymbolFilter(gname)))
      if(length(region)==0) 
        region <- reduce(genes(ensdb, filter=GeneIdFilter(gname)))
      if(length(region)==0) 
        region <- reduce(genes(ensdb, filter=TxIdFilter(gname)))
      if(length(region)==0){
        region <- reduce(genes(ensdb, filter=SymbolFilter(gname,"startsWith")))
        if(length(region)>1)
          stop("Gene not found with this exact name, and mutliple genes match ",
               "this string.")
      }
    }else{
      region <- reduce(genes(ensdb, filter=list(gene_id=gname)))
    }
    if(length(region)==0) stop("Gene/transcript not found!")
    if(length(region)>1)
      stop("Region input is ambiguous (multiple non-overlapping regions)")
    region <- strsplit(gsub("-",":",as.character(region)),":")[[1]][1:3]
  }
  stopifnot(length(region) %in% 3:4)
  region <- as.list(region)
  region[[2]] <- as.integer(region[[2]])
  region[[3]] <- as.integer(region[[3]])
  if(asGR) return(GRanges(region[[1]], IRanges(region[[2]], region[[3]])))
  region
}
