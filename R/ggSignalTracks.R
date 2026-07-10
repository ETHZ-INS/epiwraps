#' ggSignalTracks: Plot genomic signal tracks with ggplot2
#'
#' @param tracks A named list, where each element represents a track. Each track
#'   can be either 1) the path to one or multiple bigwig files (that will be 
#'   grouped), 2) an `RleList` object, or 3) a GRanges object that will be 
#'   shown as boxes.
#' @param region The region to plot, provided either as a GRanges or character.
#'   If `ensdb` is given, `region` can also be a gene name, which will be 
#'   looked up.
#' @param ensdb An optional \code{\link[ensembldb]{EnsDb}} object form which 
#' to grab transcripts. A `TxDb` object should also be supported.
#' @param extend A numeric value indicating how much to extend beyond the 
#'   `region`. If greater than 1, this indicates the number of nucleotides to 
#'   add in both directions. If smaller than or equal to 1, this indicates the
#'   proportion of the region to add on both sides. Default 0.
#' @param colors A vector of colors for each element of `tracks` (nested 
#'   elements have the same colors) which will be use to color the coverage 
#'   profiles. Colors will be recycled if necessary.
#' @param transcripts Either 'full' (full transcripts plotted), 'collapsed' 
#'   (to genes), or 'none'.
#' @param aggregation How to aggregate/show nested tracks. Either 'mean', 
#'   'heatmap' (no aggregation), or 'mean+heatmap' (both, default).
#' @param showSE Logical; whether to show the standard error on the coverage
#'   tracks of aggregated data.
#' @param nbins The number of bins in which to divide the region.
#' @param heatmap.palette A character vector specifying the colors for the 
#'   heatmap. If of length 1, is assumed to indicate the RColorBrewer palette 
#'   to use. If of length>1, the colors will be used with
#'   \code{\link[ggplot2]{scale_fill_gradientn}}.
#' @param binSummFn How to summarize date within a display bin. Either 
#'   'mean' (default) or 'max'.
#' @param sameLimits Logical; should the tracks have the same y-axis limits (and
#'   same color scale for heatmaps)?
#' @param gene_label What labels to print for genes. Either "symbol", "gene_id",
#'   "tx_name", or NULL.
#' @param trans Optional transformation of the data, either 'none' (default),
#'   'sqrt', or 'log1p'.
#' @param gene_color The genes' color.
#' @param baseTextSize The base plotting text size.
#' @param xAxis Logical; whether to plot the xAxis in the bottom panel.
#' @param coverage.linewidth Line width of the coverage plots (above ribbons).
#' @param verbose Logical; whether to print progress
#' 
#' @return A list of ggplot objects.
#'
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRanges start end seqnames restrict
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_tile geom_rect
#' @importFrom ggplot2 geom_segment geom_text scale_x_continuous margin labs
#' @importFrom ggplot2 scale_y_continuous scale_fill_distiller theme_classic 
#' @importFrom ggplot2 theme element_blank element_text unit arrow annotate
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom AnnotationFilter GRangesFilter
#' @importFrom GenomicFeatures transcripts genes
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom ensembldb genes
#' @importFrom scales comma
#' @importFrom matrixStats rowMins rowMaxs rowMedians
#' @importFrom Seqinfo seqlevels
#' @export
#' @examples
#' # we create dummy data
#' bw1 <- tempfile(fileext=".bw")
#' cov1 <- GRanges("chr1", IRanges(1L+round(c(500+rnorm(15, sd=20),
#'                                            1000*runif(20))), width=30))
#' bw2 <- tempfile(fileext=".bw")
#' cov2 <- GRanges("chr1", IRanges(1L+abs(round(c(490+rnorm(15, sd=30),
#'                                                1000*runif(20)))), width=30))
#' seqlengths(cov1) <- seqlengths(cov2) <- c("chr1"=1500)
#' rtracklayer::export.bw(coverage(cov1), bw1)
#' rtracklayer::export.bw(coverage(cov2), bw2)
#' # then we create the ggplots, and plot them:
#' pl <- ggSignalTracks(list(group=c(rep1=bw1, rep2=bw2)), region="chr1:1-1030",
#'                      aggregation="heatmap+mean")
#' patchwork::wrap_plots(pl, ncol=1, heights=c(2,1))
ggSignalTracks <- function( tracks, region, ensdb=NULL, colors="darkblue",
                            transcripts=c("full", "collapsed", "none"),
                            aggregation=c("mean+heatmap","mean","heatmap",
                                          "heatmap+mean"),
                            extend=0, showSE=FALSE, nbins=1000,
                            heatmap.palette=c("white", "blue", "black"),
                            binSummFn=c("mean", "max"), sameLimits=TRUE,
                            gene_label="symbol", trans=c("none","sqrt","log1p"),
                            gene_color="black", baseTextSize=9, xAxis=TRUE,
                            coverage.linewidth=0.2, verbose=FALSE ){
  binSummFn <- match.arg(binSummFn)
  transcripts <- match.arg(transcripts)
  aggregation <- match.arg(aggregation)
  region <- .parseRegion(region, ensdb, asGR=TRUE)
  trans <- match.arg(trans)
  stopifnot(length(extend)==1 && extend>=0)
  if(extend>0){
    if(extend<=1) extend <- extend*width(region)
    region <- .safeGRresize(region, width=2*extend+width(region), fix="center")
    isCirc <- seqinfo(region)@is_circular
    if(is.na(isCirc) || is.null(isCirc) || !isCirc){
      start(region) <- max(start(region),1L)
    }
  }
  if((2*nbins) >= width(region))
    nbins <- width(region)
  
  
  if(length(tracks)>0){
    if(is.null(names(tracks))) names(tracks) <- paste0("Track", seq_along(tracks))
  
    colors <- setNames(rep_len(colors, length(tracks)), names(tracks))
    
    if(verbose) pb <- txtProgressBar(min=0L, max=length(tracks), style=3)

    track_data_all <- lapply(seq_along(tracks), function(ti){
      gn   <- names(tracks)[ti]
      tr <- tracks[[ti]]
      if(is(tr, "GRanges"))
        return(.grangesTrack(tr, region, yname=gn, colors[ti],
                             baseTextSize=baseTextSize))
      if(inherits(tr, "RleList")){
        dat <- as.integer(Views(tr, region)[[1]][[1]])
        dat <- data.frame(pos=start(region):end(region), score=dat)
        dat <- list(.binSignal(dat, region, nbins, summFn=binSummFn))
      }else{
        bw_paths <- tracks[[ti]]
        if(is.null(names(bw_paths))) names(bw_paths) <- .getBwNames(bw_paths)
        dat <- lapply(bw_paths, function(bw){
          .binSignal(.importSingleRegionBW(bw, region), region, nbins, 
                     summFn=binSummFn)
        })
      }
      if(verbose) setTxtProgressBar(pb, ti)
      dat
    })
    if(verbose) close(pb)
    
    ymax <- NULL
    if(sameLimits){
      ymax <- unlist(lapply(track_data_all, \(x){
        if(inherits(x, "ggplot")) return(NA)
        lapply(x, \(y) max(y$score))
      }))
      if(all(is.na(ymax))){
        ymax <- NULL
      }else{
        ymax <- max(ymax, na.rm=TRUE)
      }
    }
  }
  
  panels <- list()
  
  for (gi in seq_along(tracks)){
    if(inherits(track_data_all[[gi]], "ggplot")){
      panels[[length(panels) + 1]] <- track_data_all[[gi]]
    }else{
      gn   <- names(tracks)[gi]
      dat  <- track_data_all[[gi]]
      
      if(length(dat)==1 || grepl("mean", aggregation)) {
        panels[[length(panels) + 1]] <- .coverageTrack(
          dat, gn, fill_color=colors[gi], region=region, trans=trans,
          showSE=showSE, ylim=ymax, baseTextSize=baseTextSize,
          lineWidth=coverage.linewidth)
      }
      if (length(dat)>1 && grepl("heatmap", aggregation)) {
        panels[[length(panels) + 1]] <- .heatmapTrack(
          dat, gn, palette=heatmap.palette, region=region, trans=trans,
          ymax=ymax, baseTextSize=baseTextSize)
      }
    }
  }
  
  if(!is.null(ensdb) && transcripts!="none"){
    panels$genes <- .geneTrack(ensdb, region, baseTextSize=baseTextSize,
                               collapse=transcripts=="collapsed", 
                               label=gene_label, color=gene_color, 
                               geneLabelSize=2.5)
  }
  
  if(xAxis){
    # show x-axis on last panel
    panels[[length(panels)]] <- panels[[length(panels)]] + 
      theme(axis.line.x=element_line(), axis.text.x=element_text(),
            axis.ticks.x=element_line(),
            axis.title.x=element_text()) + .bottomXLab(region)
  }
  
  panels
}


# Builds a coverage ggplot.
#
# @param binned Named list of data.frames (one per replicate), each with 
#   columns pos and score.
# @param group_name  Character label.
# @param fill_color  Colour for ribbon / line.
# @param region  The region to plot
# @param y_label y axis title
# @param ylim The upper y-axis limit (NULL for auto)
# @param trans An optional transformation
# @param baseTextSize Base plotting text size
# @param lineWidth Line width
# @param showSE  Whether to draw SE across replicates as ribbon.
.coverageTrack <- function(binned, group_name, fill_color="darkblue", region,
                           showSE=TRUE, y_label=NULL, ylim=NULL,
                           trans="none", baseTextSize=9, lineWidth=0){
  pos <- binned[[1]]$pos
  mat <- matrix(unlist(lapply(binned, \(x) x$score)), ncol=length(binned))
  n <- ncol(mat)
  mn <- rowMeans(mat, na.rm=TRUE)
  df_plot <- data.frame(pos=pos, mean=mn)
  if(n > 1){
    se <- apply(mat, 1, sd, na.rm=TRUE) / sqrt(n)
    df_plot$ymin <- mn-se
    df_plot$ymax <- mn+se
  }

  p <- ggplot2::ggplot(df_plot, aes(x=pos)) +
    ggplot2::theme_classic(base_size=baseTextSize) +
    ggplot2::theme(axis.title.x=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.line.x=element_blank(),
                   plot.margin=ggplot2::margin(2, 5, 0, 5)) +
    ggplot2::scale_x_continuous(limits=c(start(region), end(region)),
                                expand=c(0, 0)) +
    ggplot2::labs(y=ifelse(!is.null(y_label), y_label, group_name)) +
    scale_y_continuous(breaks=scales::breaks_pretty(n=2),
                       limits=c(0,ifelse(is.null(ylim), NA, ylim)),
                       transform=ifelse(trans=="none", "identity", trans))

  if(showSE && n > 1) {
    p <- p + ggplot2::geom_ribbon(aes(ymin=ymin, ymax=ymax),
                                  fill=fill_color, alpha=0.3)
      
  } else {
    p <- p +
      ggplot2::geom_ribbon(aes(ymin=0, ymax=mean),
                           fill=fill_color)
  }
  p <- p +
    ggplot2::geom_line(aes(y=mean), colour=fill_color,
                       linewidth=lineWidth) +
    ggplot2::geom_line(aes(y=0), colour=fill_color,
                     linewidth=lineWidth)
  p
}

# Build a heatmap ggplot for one group (one row per replicate).
.heatmapTrack <- function(binned, group_name,
                          palette  ="Blues",
                          region, ymax, trans="none", baseTextSize=9) {
  binned <- lapply(setNames(names(binned), names(binned)), function(nm){
    d <- binned[[nm]]
    d$sample <- factor(nm, rev(names(binned)))
    d
  })
  df_all <- do.call(rbind, binned)
  df_all$sample <- factor(df_all$sample, levels=rev(names(binned)))

  p <- ggplot(df_all, aes(x=pos, y=sample, fill=score)) + geom_tile()
  if(is(palette, "ScaleContinuous")){
    p <- p + palette
  }else if(length(palette)==1){
    p <- p + 
      scale_fill_distiller(palette=palette, direction=1, name="coverage",
                           limits=c(0, ifelse(is.null(ymax), NA, ymax)),
                           transform=ifelse(trans=="none", "identity", trans))
  }else{
    p <- p +
      scale_fill_gradientn(colours=palette, name="coverage",
                           limits=c(0, ifelse(is.null(ymax), NA, ymax)),
                           transform=ifelse(trans=="none", "identity", trans))
  }
  p <- p + ylab(group_name) +
    ggplot2::scale_x_continuous(limits=c(start(region), end(region)),
                                expand=c(0, 0)) +
    ggplot2::theme_classic(base_size=baseTextSize) +
    ggplot2::theme(axis.line=element_blank(),
                   axis.title.x =element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   #axis.title.y=element_blank(),
                   legend.key.height=ggplot2::unit(0.4, "cm"),
                   legend.key.width =ggplot2::unit(0.3, "cm"),
                   legend.text =ggplot2::element_text(size=7),
                   legend.title=ggplot2::element_text(size=7),
                   plot.margin =ggplot2::margin(0, 5, 0, 5))
  p
}

# Build a ggplot track for a GRanges object
.grangesTrack <- function(gr, yname, region, color="darkblue", baseTextSize=9){
  gr <- gr[overlapsAny(gr, region)]
  gr <- restrict(gr, start=start(region), end=end(region))
  df <- data.frame(start=start(gr), end=end(gr),
                   y=.packIntervals(start(gr), end(gr)))
  df$score <- score(gr)
  
  p <- ggplot(df, aes(xmin=start, xmax=end, ymin=y-0.4,  ymax=y+0.4))
  if(!is.null(df$score) && !all(is.na(df$score))){
    p <- p + geom_rect(aes(fill=score), colour=NA)
  }else{
    p <- p + geom_rect(fill=color, colour=NA)
  }
  p + scale_x_continuous(limits=c(start(region), end(region)), expand=c(0,0)) +
    scale_y_continuous(limits=c(0.5, max(df$y) + 0.5), breaks=NULL) +
    ylab(yname) + theme_classic(base_size=baseTextSize) + 
    theme(axis.line=element_blank(), axis.title.x=element_blank(),
          axis.text=element_blank(), axis.ticks=element_blank(),
          plot.margin=ggplot2::margin(0, 5, 0, 5))
}

# Build a ggplot gene/transcript annotation track.
#
# @param txdb      A TxDb, EnsDb, or NULL.
# @param region GRanges defining the region.
# @param xmin,xmax Numeric genomic coordinates.
# @param collapse  "gene"=one row per gene, "transcript"=one per tx.
# @param label     "symbol", "gene_id", "tx_name", or NULL.
# @param color     Arrow/exon colour.
#
# @details Requires `GenomicFeatures` for TxDb or `ensembldb` for EnsDb.
#   If neither is available or txdb is NULL, an empty placeholder is returned.
.geneTrack <- function(txdb, region, collapse=TRUE,
                       label   ="symbol",
                       color   ="#333333",
                       arrow_bins=8, baseTextSize=9, geneLabelSize=3) {
  xmax <- end(region)
  xmin <- start(region)
  
  .emptyTrack <- function(msg="") {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x=(xmin + xmax) / 2, y=0.5,
                        label=msg, size=geneLabelSize, colour="grey50") +
      ggplot2::scale_x_continuous(limits=c(xmin, xmax), expand=c(0, 0)) +
      ggplot2::scale_y_continuous(limits=c(0, 1)) +
      ggplot2::theme_classic(base_size=baseTextSize) +
      ggplot2::theme(axis.title.y =element_blank(),
                     axis.text.y  =element_blank(),
                     axis.ticks.y =element_blank(),
                     plot.margin  =ggplot2::margin(0, 5, 2, 5))
  }

  if (is.null(txdb)) return(.emptyTrack("No gene annotation provided"))

  is_ensdb <- inherits(txdb, "EnsDb")
  if (is_ensdb) {
    filter  <- AnnotationFilter::GRangesFilter(region)
    exons   <- exons(txdb, filter=filter, 
                     columns=c("gene_id", "gene_name", "tx_id"))
    if(!collapse){
      txs <- transcripts(txdb, filter=filter,
                                      columns=c("gene_id", "gene_name",
                                                  "tx_id"))
    }else{
      exons$tx_id <- exons$gene_name
      txs <- genes(txdb, filter=filter, columns=c("gene_id", "gene_name"))
    }
  } else {
    exons <- exons(txdb, columns=c("gene_id", "tx_name"),
                   filter=list("gene_chrom"=
                                 as.character(GenomicRanges::seqnames(region))))
    exons <- IRanges::subsetByOverlaps(exons, region)
    if(!collapse){
      exons$tx_id <- exons$tx_name
      txs   <- GenomicFeatures::transcripts(txdb,
                  columns =c("gene_id", "tx_name"),
                  filter  =list(
                    "tx_chrom"=as.character(
                      GenomicRanges::seqnames(region))))
      txs   <- IRanges::subsetByOverlaps(txs, region)
      txs$tx_id <- txs$tx_name
    }else{
      exons$tx_id <- sapply(exons$gene_id, `[`, 1)
      txs <- GenomicFeatures::genes(txdb,
               filter=list("tx_chrom"=as.character(seqnames(region))))
    }
  }
  
  if(collapse){
    tmp <- reduce(GRangesList(split(exons, exons$gene_name)))
    exons <- unlist(tmp)
    exons$tx_id <- rep(names(tmp), lengths(tmp))
  }

  if (length(txs) == 0) return(.emptyTrack("No gene in region"))

  tx_df <- as.data.frame(txs, row.names=NULL)
  tx_df$gname <- if (is_ensdb) tx_df$gene_name else
    sapply(tx_df$gene_id, `[`, 1)
  if(collapse){
    tx_df$tx_id <- tx_df$gname
  }else{
    tx_df$tx_id <- if (is_ensdb) tx_df$tx_id else tx_df$tx_name
  }
  
  # VERY MUCH CLAUDE CODE FROM HERE ON...

  # Assign a y-level (pack transcripts so they don't overlap)
  tx_df <- tx_df[order(tx_df$start), ]
  tx_df$ylevel <- .packIntervals(tx_df$start, tx_df$end)
  row.names(tx_df) <- NULL

  # Exon data.frame - convert with row.names=NULL to avoid duplicate issues
  ex_df <- as.data.frame(exons, row.names=NULL)
  ex_df$tx_id <- if (is_ensdb) ex_df$tx_id else ex_df$tx_name
  # Clip to region
  ex_df$start <- pmax(ex_df$start, start(region))
  ex_df$end   <- pmin(ex_df$end,   end(region))
  tx_df$start_cl <- pmax(tx_df$start, start(region))
  tx_df$end_cl   <- pmin(tx_df$end,   end(region))
  ex_df$exon_id <- NULL

  # Merge ylevel into exons
  ex_df <- merge(ex_df, tx_df[, c("tx_id", "ylevel")], by="tx_id",
                 all.x=TRUE)

  # direction arrows along introns
  arrow_df <- lapply(seq_len(nrow(tx_df)), function(i) {
    tx  <- tx_df[i, ]
    xs  <- seq(tx$start_cl, tx$end_cl, length.out=arrow_bins + 1)
    xs  <- xs[-1]  # midpoints of bins
    if (tx$strand == "-") xs <- rev(xs)
    data.frame(x=head(xs, -1), xend=tail(xs, -1),
               y=tx$ylevel, yend=tx$ylevel,
               tx_id=tx$tx_id)
  })
  arrow_df <- do.call(rbind, arrow_df)

  # label data.frame
  lbl_df <- tx_df
  lbl_df$lbl <- switch(label,
    symbol = if("gname"  %in% names(lbl_df)) lbl_df$gname  else lbl_df$tx_id,
    gene_id = if("gene_id" %in% names(lbl_df)){
      sapply(lbl_df$gene_id, `[`, 1) }else{ lbl_df$tx_id },
    tx_name = lbl_df$tx_id,
    NULL
  )
  lbl_df$lbl_x <- (lbl_df$start_cl + lbl_df$end_cl) / 2

  n_levels <- max(tx_df$ylevel, na.rm=TRUE)

  p <- ggplot2::ggplot() +
    # thin backbone line per transcript
    ggplot2::geom_segment(data=tx_df,
      aes(x=start_cl, xend=end_cl,
                   y=ylevel, yend=ylevel),
      colour=color, linewidth=0.4) +
    # direction arrows
    ggplot2::geom_segment(data=arrow_df,
      aes(x=x, xend=xend, y=y, yend=yend),
      colour=color, linewidth=0.3,
      arrow=ggplot2::arrow(length=ggplot2::unit(0.08, "cm"),
                             type  ="open")) +
    # exon boxes
    ggplot2::geom_rect(data=ex_df,
      aes(xmin=start, xmax=end,
                   ymin=ylevel - 0.25, ymax=ylevel + 0.25),
      fill=color, colour=NA) +
    ggplot2::scale_x_continuous(limits=c(start(region), end(region)),
                                expand=c(0, 0),
                                labels=scales::comma) +
    ggplot2::scale_y_continuous(limits=c(0.5, n_levels + 0.5),
                                breaks =NULL) +
    ggplot2::theme_classic(base_size=baseTextSize) +
    ggplot2::theme(axis.line=element_blank(),
                   axis.title.x =element_blank(),
                   axis.text=element_blank(),
                   axis.ticks=element_blank(),
                   plot.margin =ggplot2::margin(0, 5, 2, 5))

  if(is_ensdb) p <- p + ylab(paste0("ensembl\n", ensemblVersion(txdb)))
  
  if (!is.null(label)) {
    p <- p +
      ggplot2::geom_text(data=lbl_df,
        aes(x=lbl_x, y=ylevel + 0.35, label=lbl),
        size=geneLabelSize, colour=color, fontface="italic",
        hjust=0.5, vjust=0)
  }
  p
}

.bottomXLab <- function(region){
  xlab(paste0(as.character(seqnames(region)),
                " : ", scales::comma(start(region)), " - ", 
                scales::comma(end(region))))
}

# Pack intervals onto as few non-overlapping levels as possible.
# Returns an integer vector of y-levels.
# This is Claude Code... thanks to whoever coded something similar!
.packIntervals <- function(starts, ends){
  n <- length(starts)
  levels <- integer(n)
  ends_of_levels <- c()
  for (i in seq_len(n)) {
    placed <- FALSE
    for (lv in seq_along(ends_of_levels)) {
      if (starts[i] > ends_of_levels[lv] + 1) {
        levels[i] <- lv
        ends_of_levels[lv] <- ends[i]
        placed <- TRUE
        break
      }
    }
    if(!placed){
      ends_of_levels <- c(ends_of_levels, ends[i])
      levels[i] <- length(ends_of_levels)
    }
  }
  levels
}


.getBwNames <- function(paths, ext.regex="\\.bw$|\\.bigwig$"){
  if(!is.null(ext.regex))
    paths <- gsub(ext.regex, "", paths, ignore.case=TRUE)
  bn <- basename(paths)
  if(!any(duplicated(bn))) return(bn)
  if(length(unique(bn))==1 && !all(paths==dirname(paths))){
    return(.getBwNames(dirname(paths), NULL))
  }
  paste0("rep", seq_along(paths))
}

.importSingleRegionBW <- function(bw_path, region) {
  gr  <- rtracklayer::import.bw(bw_path, which=region)
  if(length(gr) == 0) return(data.frame(pos=integer(0), score=numeric(0)))
  gr <- GenomicRanges::restrict(gr, start=GenomicRanges::start(region),
                                end=GenomicRanges::end(region))
  rows   <- lapply(seq_along(gr), function(i)
    data.frame(pos=start(gr)[i]:end(gr)[i], score=gr$score[i]))
  do.call(rbind, rows)
}

# Bin a per-base data.frame to `nbins` bins within the region.
.binSignal <- function(df, region, nbins=1000L, summFn=c("mean","max")){
  if(nrow(df)==0)
    return(data.frame(pos=seq(start(region), end(region), length.out=nbins),
                      score=0))
  breaks  <- seq(start(region), end(region), length.out=nbins+1)
  mids    <- (breaks[-1] + breaks[-(nbins+1)])/2
  df$bin  <- cut(df$pos, breaks=breaks, include.lowest=TRUE, labels=FALSE)
  summFn <- switch(match.arg(summFn), max=max, mean=mean)
  agg <- tapply(df$score, df$bin, summFn, na.rm=TRUE)
  out <- data.frame(pos=mids, score=0)
  out$score[as.integer(names(agg))] <- agg
  out
}
