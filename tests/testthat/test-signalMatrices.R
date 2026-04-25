bw <- system.file("extdata/example_atac.bw", package="epiwraps")
regions <- system.file("extdata/example_peaks.bed", package="epiwraps")
regions <- importBedlike(regions)

checkESE <- function(m, enc=1L, enr=NULL, wAssay=1L){
  expect_true(is(m, "EnrichmentSE") && ncol(m)==enc)
  if(!is.null(enr)) expect_true(nrow(m)==enr)
  expect_all_true(unlist(lapply(assay(m,wAssay),
                                \(x) is.finite(as.numeric(x)))))
}

test_that("signal2Matrix works with bw files", {
  m <- signal2Matrix(bw, regions, scaling=2, w=20)
  checkESE(m, enr=length(regions))
  m <- signal2Matrix(bw, regions, type = "scaled", smooth=TRUE)
  checkESE(m, enr=length(regions))
  h <- plotEnrichedHeatmaps(m)
})

m_ref <- signal2Matrix(list(test=bw), regions)

test_that("signal2Matrix works with RleList", {
  rle <- rtracklayer::import(bw, as="RleList")
  m <- signal2Matrix(list(test=rle), regions)
  expect_true(identical(m,m_ref))
})

test_that("signal2Matrix works with GRanges", {
  gr <- rtracklayer::import(bw)
  m <- signal2Matrix(list(test=gr), regions)
  cc <- cor(as.numeric(assay(m)[,1]),as.numeric(assay(m_ref)[,1]))
  expect_true(cc>0.99)
})

bam <- system.file("extdata", "ex1.bam", package="Rsamtools")

test_that("signal2Matrix works with Bam", {
  m <- signal2Matrix(bam, as(c("seq1:500-800"), "GRanges"))
  checkESE(m, enr = 1)
})


data("exampleESE")

test_that("Clustering works", {
  rowData(exampleESE)$cluster <<- clusterSignalMatrices(exampleESE, k=3,
                                                        scaleRows=TRUE)
  expect_equal(length(table(rowData(exampleESE)$cluster)), 3L)
})

test_that("plotEnrichedHeatmaps works", {
  pdf(NULL)
  expect_no_error(draw(plotEnrichedHeatmaps(exampleESE, multiScale=TRUE,
                        scale_rows="global", minRowVal=1, row_split="cluster",
                        mean_color=c("1"="red", "2"="blue", "3"="black"),
                        colors=list("darkblue","darkgreen","darkred"))))
  dev.off()
})


test_that("Melting works", {
  d <- meltSignals(exampleESE, splitBy = "cluster")
  expect_equal(nrow(d), 720L)
  expect_all_true(unlist(lapply(d, \(x) sum(is.na(x))==0)))
  d <- meltSignals(m_ref)
})



test_that("Merging works", {
  test <- mergeSignalMatrices(exampleESE[,1:2], "median")
  expect_all_true(dim(test)==c(150L,80L))
  expect_all_true(is.finite(as.numeric(test)))
  test <- mergeSignalMatrices(exampleESE[,1:2], "mean")
  expect_all_true(dim(test)==c(150L,80L))
  ese2 <- ml2ESE(list(merged=test), rowRanges=rowRanges(exampleESE))
  checkESE(ese2, enr = nrow(exampleESE))
  ese2 <- cbind(exampleESE, ese2)
  checkESE(ese2, 4L)
})


test_that("ESE methods works", {
  expect_no_error(show(exampleESE))
  expect_no_error(showTrackInfo(exampleESE))
  ml <- getSignalMatrices(exampleESE)
  expect_true(is.list(ml) && length(ml)==3)
  expect_all_true(is.finite(as.numeric(score(exampleESE))))
})

test_that("Normalization of signalMatrices works", {
  test <- renormalizeSignalMatrices(exampleESE, method="border")
  checkESE(test, 3L)
  test <- renormalizeSignalMatrices(exampleESE)  
  checkESE(test, 3L)
  test <- renormalizeSignalMatrices(exampleESE, method="manual",
                                    scaleFactors=c(1,1.2,0.8))
  checkESE(test, 3L)
})


test_that("getNormFactors works with bw", {
  nf <- suppressWarnings(getNormFactors(c(bw,bw), useSeqLevels="1"))
  expect_all_true(is.finite(nf))
})

test_that("getNormFactors works with bam", {
  nf <- suppressWarnings(getNormFactors(c(bam,bam), paired=FALSE, nwind = 10L))
  expect_all_true(is.finite(nf))
})

test_that("matrix resizing works", {
  out <- resizeMatrix(assay(exampleESE)[,1], ndim = c(30,10))
  expect_all_true(is.finite(as.numeric(out)))
  expect_all_true(dim(out)==c(30,10))
})
