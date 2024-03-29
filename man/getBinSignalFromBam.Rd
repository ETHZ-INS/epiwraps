% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signal2Matrix.R
\name{getBinSignalFromBam}
\alias{getBinSignalFromBam}
\title{getBinSignalFromBam}
\usage{
getBinSignalFromBam(
  filepath,
  regions,
  cuts = FALSE,
  RPM = TRUE,
  paired = TRUE,
  ...,
  flgs = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE)
)
}
\arguments{
\item{filepath}{The path to the (indexed) bam file}

\item{regions}{A `GRanges` of the regions/positions around which to plot}

\item{cuts}{Whether to count cuts (e.g. beginning/end of fragments) rather
than coverage (ignored unless the input are bam files)}

\item{RPM}{Whether to perform RPM normalization (for bam input)}

\item{paired}{Whether to consider whole fragments}

\item{flgs}{Flags for bam reading}
}
\value{
A matrix
}
\description{
This is a wrapper around \code{\link[genomation]{ScoreMatrix}} to enable BAM
support in \code{\link{signal2Matrix}}.
}
