% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization.R
\name{renormalizeSignalMatrices}
\alias{renormalizeSignalMatrices}
\title{renormalizeSignalMatrices}
\usage{
renormalizeSignalMatrices(
  ml,
  method = c("border", "top", "manual"),
  trim = NULL,
  fromAssay = "input",
  toAssay = NULL,
  nWindows = NULL,
  scaleFactors = NULL,
  ...
)
}
\arguments{
\item{ml}{A named matrix list or EnrichmentSE object as produced by 
\code{\link{signal2Matrix}}.}

\item{method}{Either "border" or "top" (see details below).}

\item{trim}{Quantiles trimmed at each extreme before calculating 
normalization factors.}

\item{fromAssay}{Assay to use (ignored unless `ml` is an EnrichmentSE 
object), defaults to the first assay.}

\item{toAssay}{Assay in which to store the normalized data (ignored unless 
`ml` is an EnrichmentSE object). By default an assay name will be set based
on the normalization method used.}

\item{scaleFactors}{A numeric vector of same length as `ml`, 
indicating the scaling factors by which to multiply each matrix.
Alternatively, a numeric matrix with a number of rows equal to the length 
of `ml`, and two columns indicating the alpha and beta arguments of a 
s3norm normalization. Ignored unless `method="manual"`.}
}
\value{
Either a renormalized list of signal matrices or, if `ml` was an 
  `EnrichmentSE` object, the same object with an additional normalized
  assay automatically put at the front.
}
\description{
Renormalizes a list of signal matrices or an EnrichmentSE object.
}
\details{
* `method="border"` works on the assumption that the left/right borders of the 
matrices represent background signal which should be equal across samples. As
a result, it will work only if 1) the left/right borders of the matrices are 
sufficiently far from the signal (e.g. peaks) to be chiefly noise, and 
2) the signal-to-noise ratio is comparable across tracks/samples.
* `method="top"` instead works on the assumption that the highest signal should
be the same across tracks/samples.
By default, extreme values are trimmed before establishing either kind of 
normalization factor. The proportion trimmed can be set using the `trim` 
argument, and is by default 10% of the non-zero values.
* `method="manual"` enables the use of independently computed normalization
factors, for instance obtained through \code{\link{getNormFactors}}.
}