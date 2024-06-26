% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization.R
\name{getNormFactors}
\alias{getNormFactors}
\alias{bwNormFactors}
\title{getNormFactors : estimate normalization factors from genomic signal files}
\usage{
getNormFactors(
  x,
  method = c("background", "SES", "enriched", "top", "MAnorm", "S3norm", "2cLinear"),
  wsize = 10L,
  nwind = 20000L,
  peaks = NULL,
  trim = 0.02,
  useSeqLevels = NULL,
  paired = NULL,
  ...,
  verbose = TRUE
)

bwNormFactors(x, ...)
}
\arguments{
\item{x}{A vector of paths to bigwig files, to bam files, or alternatively a 
list of coverages in RleList format. (Mixed formats are not supported)}

\item{method}{Normalization method (see details below).}

\item{wsize}{The size of the random windows. If any of the bigwig files 
records point events (e.g. insertions) at high resolution (e.g. nucleotide),
use a lot (e.g. >10k) of small windows (e.g. `wsize=10`), as per default
settings. Otherwise the process can be lightened by using fewer bigger 
windows.}

\item{nwind}{The number of random windows}

\item{peaks}{A list of peaks (GRanges) for each experiment in `x`, or a
vector of paths to such files, or a single GRanges of unified peaks to use 
(e.g. for top/MAnorm).}

\item{trim}{Amount of trimming when calculating means.}
}
\value{
A vector of normalization factors, or for the 'S3norm' and '2cLinear'
  methods, a numeric matrix with a number of rows equal to the length 
  of `x`, and two columns indicating the alpha and beta terms.

A vector of normalization factors or, for methods 'S3norm' and 
  '2cLinear', a matrix of per-sample normalization parameters.
}
\description{
Estimates normalization factors for a set of samples (i.e. bam/bigwig files).
}
\details{
The function computes per-sample (bigwig or bam file) normalization factors
using one of the following methods:
* The 'background' or 'SES' normalization method (they are synonyms here)
(Diaz et al., Stat Appl Gen et Mol Biol, 2012) assumes that the background
noise should on average be the same across experiments, an assumption that 
works well in practice when there are not very large differences in 
signal-to-noise ratio. The implementation uses the trimmed mean number of 
reads in random windows with non-zero counts.
* The 'MAnorm' approach (Shao et al., Genome Biology 2012) assumes that 
regions that are commonly enriched (i.e. common peaks) in two experiments 
should on average have the same signal in the two experiments. If the data is
strictly positive (e.g. counts or CPMs), TMM normalization (Robinson et al.,
Genome Biology 2010) is then employed on the common peaks. When done on more 
than two samples, the sample with the highest number of overlaps with other
samples is used as a reference.
* The 'enriched' approach assumes that enriched regions are on average 
similarly enriched across samples. Contrarily to 'MAnorm', these regions do 
not need to be in common across samples/experiments. Note that trimming via
the `trim` parameter is performed before calculating means.
* The 'top' approach assumes that the maximum enrichment (after trimming 
according to the `trim` parameter) in peaks is the same across 
samples/experiments.
* The 'S3norm' (Xiang et al., NAR 2020) and '2cLinear' methods try to 
normalize both enrichment and background simultaneously. S3norm does this in 
a log-linear fashion (as in the publication), while '2cLinear' does it on 
the original scale.

Several of the methods rely on random regions to estimate background levels.
This means that they will not be entirely deterministic, although normally 
quite reproducible with experiments that are not very sparse. For fully 
reproducible results, be sure to set the random seed.

When the data is very sparse (e.g. low sequencing depth, or compiling single
nucleotide hits rather than fragment coverage), it might be necessary to 
increase the `wsize` and `nwind` to get more robust estimates.
}
\section{Functions}{
\itemize{
\item \code{bwNormFactors()}: deprecated in favor of getNormFactors

}}
\examples{
# we get an example bigwig file, and use it twice:
bw <- system.file("extdata/example_atac.bw", package="epiwraps")
bw <- c(sample1=bw, sample2=bw)
# we indicate to use only chr1, because it's the only one in the example file
getNormFactors(bw, useSeqLevels="1")
}
