% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc.R
\name{plotCovStats}
\alias{plotCovStats}
\title{plotCovStats}
\usage{
plotCovStats(qc, labels = "AUTO", show.legend = TRUE)
}
\arguments{
\item{qc}{A list of coverage statistics, as produced by 
\code{\link{getCovStats}}.}

\item{labels}{Passed to \code{\link[cowplot]{plot_grid}}.}
}
\value{
A grid object to be plotted.
}
\description{
Plots coverage statistics, such as as fingerprint plot.
}
\examples{
# we use an example bigwig file
bwf <- system.file("extdata/example_atac.bw", package="epiwraps")
# because most of the file is empty, we'll exclude some of the ranges
cs <- getCovStats(bwf, exclude=GRanges("1", IRanges(1, 4300000)))
plotCovStats(cs)
}
