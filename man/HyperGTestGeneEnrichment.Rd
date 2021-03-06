% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cAMARETTO_Results.R
\name{HyperGTestGeneEnrichment}
\alias{HyperGTestGeneEnrichment}
\title{GmtFromModules}
\usage{
HyperGTestGeneEnrichment(
  gmtfile1,
  gmtfile2,
  runname1,
  runname2,
  NrCores,
  ref.numb.genes = 45956
)
}
\arguments{
\item{gmtfile1}{A gmtfilename that you want to compare}

\item{gmtfile2}{A second gmtfile to compare with.}

\item{runname1}{name of the first dataset.}

\item{runname2}{name of the second dataset.}

\item{NrCores}{Number of cores for parallel computing.}

\item{ref.numb.genes}{The reference number of genes.}
}
\value{
Creates resultfile with p-values and padj when
comparing two gmt files with a hyper geometric test.
}
\description{
GmtFromModules
}
\examples{
try(
HyperGTestGeneEnrichment(gmtfile1, gmtfile2,
runname1, runname2,NrCores, ref.numb.genes=45956)
)
}
