bam <- system.file("extdata", "ex1.bam", package="Rsamtools")

test_that("peak calling from bam works",{
  peaks <- callPeaks(bam, paired=TRUE, nChunks=1)
  expect_true(is(peaks, "GRanges"))
  peaks <- callPeaks(bam, paired=FALSE, nChunks=1)
  expect_true(is(peaks, "GRanges"))
})

test_that("peak calling from cov works",{
  bwf <- system.file("extdata/example_rna.bw", package="epiwraps")
  rle <- rtracklayer::import(bwf, as="RleList")
  peaks <- callPeaks(rle, paired=FALSE)
  expect_true(is(peaks, "GRanges"))
})

test_that("peakCountsFromBAM works",{
  peaks <- GRanges(c("seq1","seq1","seq2"), IRanges(c(400,900,500), width=100))
  se <- peakCountsFromBAM(bam, peaks, paired=FALSE, verbose=FALSE)
  expect_all_true(dim(se)==c(3,1))
  expect_all_true(as.numeric(assay(se))==c(145L, 144L, 180L))
  se <- peakCountsFromBAM(bam, peaks, paired=TRUE, verbose=FALSE)
  expect_all_true(as.numeric(assay(se))==c(168,164,207))
})


frags <- tempfile(fileext = ".tsv")
d <- data.frame(chr=rep(letters[1:2], each=10), start=rep(100*(1:10),2))
d$end <- d$start + 15L
set.seed(123)
d$cell <- paste0("barcode",sample.int(3, nrow(d), replace=TRUE))
write.table(d, frags, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
frags <- Rsamtools::bgzip(frags)
Rsamtools::indexTabix(frags, format = "bed")
regions <- GRanges(c("a","b"), IRanges(400,width=300))

test_that("peakCountsFromFrags works",{
  se <- peakCountsFromFrags(frags, regions)
  expect_all_true(dim(se)==c(2,3))
  expect_all_true(as.integer(assay(se))==c(0L, 1L, 2L, 1L, 1L, 1L))
  bcmap <- setNames(c("PB1","PB1","PB2"),paste0("barcode",1:3))
  se <- peakCountsFromFrags(frags, regions, barcodes=bcmap)
  expect_all_true(dim(se)==c(2,2))
})

test_that("reduceWithResplit works",{
  gr <- GRanges("1", IRanges(c(100,120,140,390,410,430,120),
                             width=rep(c(200,520),c(6,1))))
  redGr <- reduceWithResplit(gr, softMaxSize=100)
  expect_equal(length(redGr), 2L)
})

set.seed(123)
grl <- lapply(c(A=10,B=20,C=30), FUN=function(x){
  gr <- GRanges("seq1", IRanges(runif(x,1,1000), width=20))
  seqlengths(gr)["seq1"] <- 1200
  gr
})

test_that("overlap functions work", {
  o <- regionOverlaps(grl, mode="pairwise", returnValues=TRUE)
  expect_true(all(is.finite(o)))
  o <- regionOverlaps(grl, returnValues=TRUE)
  expect_true(all(is.finite(o)))
  o <- regionOverlaps(grl, returnValues=TRUE, colorBy="jaccard")
  expect_true(all(is.finite(o)))
  o <- regionsToUpset(grl)
  expect_all_true(colSums(o)==c(8,13,14))
  o <- regionsToUpset(grl, reference = "disjoin")
  expect_all_true(colSums(o)==c(30,53,69))
})


test_that("plotSignalTracks works", {
  bw <- system.file("extdata/example_atac.bw", package="epiwraps")
  pdf(NULL)
  expect_no_error(plotSignalTracks(list(track1=bw),
                                   region="8:22165140-22212326"))
  dev.off()
})
