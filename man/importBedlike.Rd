% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{importBedlike}
\alias{importBedlike}
\title{importBedlike}
\usage{
importBedlike(x, ...)
}
\arguments{
\item{x}{The path to a bed or bed-like file (can be 
gzipped)}

\item{...}{passed to \code{\link[data.table]{fread}}}
}
\value{
A `GRanges` object
}
\description{
Imports a bed-like file as a GRanges object. Uses
`rtracklayer` import functions if possible, and falls
back onto an import that's not format-committed 
otherwise.
}
\examples{
# example bed file:
filepath <- system.file("extdata/example_peaks.bed", 
                        package="epiwraps")
b <- importBedlike(filepath)
}
