% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ESE.R
\name{addAssayToESE}
\alias{addAssayToESE}
\title{addAssayToESE}
\usage{
addAssayToESE(x, a, name = "normalized", replace = TRUE)
}
\arguments{
\item{x}{An object of class `EnrichmentSE`, as produced by 
\code{\link{signal2Matrix}}.}

\item{a}{The assay to add, e.g. a list of normalizedMatrix objects}

\item{name}{}

\item{replace}{Logical, whether to replace any existing assay of the same 
name (default TRUE). If FALSE and the assay already existed, the new assay 
name is given a suffix.}
}
\value{
`x` with the added/updated assay.
}
\description{
Adds an assay of signal matrices to an existing `EnrichmentSE` object.
}
\examples{
# we first get an EnrichmentSE object:
data(exampleESE)
# then we will create a new assay which is simply sqrt-transformed, and add 
# it back in the object
newAssay <- lapply(getSignalMatrices(x), sqrt)
exampleESE <- addAssayToESE(exampleESE, newAssay, named="sqrt")
}
