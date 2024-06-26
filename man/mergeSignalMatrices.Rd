% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signalMatrices.R
\name{mergeSignalMatrices}
\alias{mergeSignalMatrices}
\title{mergeSignalMatrices: aggregates two or more signal matrices.}
\usage{
mergeSignalMatrices(ml, aggregation = c("mean", "sum", "median"), assay = 1L)
}
\arguments{
\item{ml}{A named list of signal matrices or an EnrichmentSE object as 
produced by \code{\link{signal2Matrix}}}

\item{aggregation}{The method to aggregate matrices}

\item{assay}{Assay to use (ignored unless `ml` is an ESE object), defaults to
the first assay.}
}
\value{
A single `normalizedMatrix` object.
}
\description{
mergeSignalMatrices: aggregates two or more signal matrices.
}
\examples{
# we first get an EnrichmentSE object:
data(exampleESE)
# we merge the two tracks:
merged <- mergeSignalMatrices(exampleESE)
# we could then plot the merge (not run):
# plotEnrichedHeatmaps(merged)
}
