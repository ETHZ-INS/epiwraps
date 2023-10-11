#' refragment
#'
#' Uses the cut profiles to re-fragment in silico fragments that contain more 
#'   than a nucleosome.
#'
#' @param bam Path to a paired-end bam file, or a `GRanges` of fragments.
#' @param minSize Minimum fragment size
#' @param nfr Maximum size of nucleosome-free fragments
#' @param nuc Range (min & max) of mono-nucleosome fragment sizes
#' @param binSize Bin size for cut profile
#' @param verbose Logical; whether to print progress messages
#' @param BPPARAM `BiocParallel` param for multithreading.
#'
#' @return A list of processed fragment `GRanges`, including nucleosome-free,
#'   mono-nucleosome, and ambiguous fragments.
#' @export
refragment <- function(bam, minSize=20L, nfr=120L, nuc=c(145,190), binSize=10L,
                       verbose=TRUE, BPPARAM=BiocParallel::SerialParam()){
  if(is(bam, "GRanges")){
    a <- bam
  }else{
    if(verbose) message("Reading alignments")
    a <- GRanges(readGAlignmentPairs(bam))
  }
  nFrags <- length(a)
  if(verbose) message("Computing breaks coverages")
  co <- tileRle(coverage(.align2cuts(a)), bs=as.integer(binSize))
  
  wi <- width(a)
  
  # initial split of fragments
  fout <- list(
    nf=a[wi<=nfr & wi>minSize],
    nuc=a[wi>=nuc[1] & wi<=nuc[2]],
    other=a[c()]
  )
  fam <- a[(wi>nfr & wi<nuc[1]) | wi>nuc[2]]
  n <- c(sum(wi<minSize), lengths(fout)[1:2], length(fam))
  pc <- round(100*n/length(a))
  if(verbose) message("Initial composition:
    below min frag size: ", n[1], " (",pc[1],"%)
    nucleosome-free: ", n[2], " (",pc[2],"%)
    mono-nucleosome: ", n[3], " (",pc[3],"%)
    ambiguous: ", n[4], " (",pc[4],"%) --> will re-fragment
    ")
  rm(a)
  gc(verbose=FALSE)
  
  wn <- width(fout$nuc)
  # Force single thread when dealing with few fragments:
  if(length(fam)<1000) BPPARAM <- SerialParam()
  # Split work by chrs
  fam <- split(ranges(fam), seqnames(fam))
  names(chrs) <- chrs <- names(fam)[lengths(fam)>0]
  fr <- bplapply(chrs, BPPARAM=BPPARAM,
                 FUN=function(x){
              .refragment(fam[[x]], co[[x]], wnuc=wn, minSize=minSize, nfr=nfr)
                 })
  fout$nf <- c(fout$nf, .viewl2gr(lapply(fr, FUN=function(x) x[["frag.nf"]])))
  fout$nuc <- c(fout$nuc, .viewl2gr(lapply(fr, FUN=function(x) x[["frag.nuc"]])))
  fout$other <- c(fout$other, .viewl2gr(lapply(fr, FUN=function(x) x[["unfrag"]])))

  n <- lengths(fout)[c("nf", "nuc")]
  if(verbose) message("Counts after re-fragmentation:
    nucleosome-free: ", n[1], "
    mono-nucleosome: ", n[2], "
    ")
  
  lapply(fout, sort)
}

# called by `refragment` above
.refragment <- function(fr, co, wnuc=wn, minSize=20L, nfr=120L, thres=10L){
  stopifnot(is(fr, "IRanges"))
  stopifnot(is(co, "Rle"))
  nuc <- range(wnuc)
  ## cut profile along fragments (ignoring first&last nt)
  v <- Views(co, resize(fr, width=width(fr)-2L, fix="center"))
  ## keep only frags with overlaping cuts
  wHasCut <- which(max(v)>0L)
  v <- RleList(v[wHasCut])
  ## find maxima
  vmax <- max(v)
  pos <- which(v==vmax)
  ## lengths on both sides of the maxima
  len1 <- width(fr)[wHasCut]-pos
  len2 <- pos+1L
  
  pf <- function(x){ # rough probability of fragments
    x <- unlist(x)
    x <- as.integer(x>=minSize) * # frags below min size not allowed
      matrixStats::rowMaxs(cbind( # frags should be either:
        1-abs(pnorm(x, mean(wnuc), sd(wnuc))-0.5)*2,  # within mononuc sizes, or
        as.numeric(x > (minSize+nuc[2]))/2, # above it (further cuts), or
        as.numeric(x < nfr)/3 )) # within nuc-free sizes
    sqrt(x)
  }
  ## combine prob of both fragments with breakpoint counts
  p <- relist(as.integer(round(100*pf(len1)*pf(len2)*
                                 rep(ppois(vmax,2.5),lengths(pos)))), pos)
  
  ## get max breakpoint count
  p.max <- max(p)
  p.w <- which.max(p)
  
  ## prepare return:
  names(otypes) <- otypes <- c("unfrag","frag.amb","frag.nuc","frag.nf")
  ret <- lapply(otypes, FUN=function(x) IRanges())
  
  if(length(wB <- which(p.max>thres))==0){
    # nothing to refragment
    ret$unfrag <- fr
    return(ret)
  }
  ## indices of frags to be refragmented, rel. to original `fr`:
  toFrag <- wHasCut[wB]
  ret$unfrag <- fr[!toFrag]
  fr <- fr[toFrag]
  ## recover the position in the read of the selected maximum:
  p <- p.w[wB]
  p <- unlist(pos[wB])[p+head(c(0L,cumsum(lengths(pos[wB]))),length(p))]
  ## create the sub-fragments:
  fr <- c(IRanges(start(fr), start(fr)+p), IRanges(start(fr)+p+1L, end(fr)))
  ## classify and return:
  w <- width(fr)
  wNuc <- w>=nuc[1] & w<=nuc[2]
  wNf <- w<=nfr
  ret$frag.nuc <- fr[which(wNuc)]
  ret$frag.nf <- fr[which(wNf)]
  frag.amb <- fr[which(!wNuc & !wNf)]
  # loop until no newly-refragmented frag is ambiguous
  while(length(frag.amb)>0){
    o <- .refragment(frag.amb, co, wnuc=wnuc, minSize=minSize, nfr=nfr, 
                     thres=thres)
    for(f in c("unfrag","frag.nuc","frag.nf")) ret[[f]] <- c(ret[[f]], o[[f]])
    frag.amb <- o$frag.amb
  }
  lapply(ret, sort)
}
