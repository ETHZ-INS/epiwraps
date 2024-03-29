% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overlaps.R
\name{regionUpset}
\alias{regionUpset}
\title{regionUpset}
\usage{
regionUpset(
  x,
  reference = c("reduce", "disjoin"),
  returnList = FALSE,
  ignore.strand = FALSE,
  maxgap = -1L,
  minoverlap = 0L,
  ...
)
}
\arguments{
\item{x}{A named list of genomic ranges (or paths to bed files)}

\item{reference}{The method for creating the reference windows ('reduce' or
'disjoin'). Alternatively, a `GRanges` object of reference windows.}

\item{returnList}{Logical; whether to return the list instead of plotting.}

\item{ignore.strand}{Logical; whether to ignore strands when computing 
overlaps (default FALSE). Strand information is ignored if either of the 
compared sets of regions is unstranded.}

\item{maxgap}{Maximum gap between regions to count as an overlap (see 
\code{\link[GenomicRanges]{findOverlaps-methods}}).}

\item{minoverlap}{Minimum overlap to count as a match (see 
\code{\link[GenomicRanges]{findOverlaps-methods}}).}

\item{...}{Further plotting arguments passed to \code{\link[UpSetR]{upset}}.}
}
\value{
A plot
}
\description{
A wrapper around \code{\link[UpSetR]{upset}} for comparing multiple sets of
genomic ranges.
}
\examples{
# random list of GRanges:
grl <- lapply(c(A=10,B=20,C=30), FUN=function(x){
  GRanges("seq1", IRanges(runif(x,1,1000), width=20))
})
regionUpset(grl)
}
