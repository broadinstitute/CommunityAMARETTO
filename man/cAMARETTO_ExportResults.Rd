% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cAMARETTO_DataReadDownload.R
\name{cAMARETTO_ExportResults}
\alias{cAMARETTO_ExportResults}
\title{cAMARETTO_ExportResults}
\usage{
cAMARETTO_ExportResults(
  cAMARETTOresults,
  cAMARETTOnetworkM,
  cAMARETTOnetworkC,
  output_address = "./"
)
}
\arguments{
\item{cAMARETTOresults}{The output of the Results function.}

\item{cAMARETTOnetworkM}{The output of the Module Network function.}

\item{cAMARETTOnetworkC}{The output of the Identify Communities function.}

\item{output_address}{The address where the zipfile will be generated in.}
}
\value{
a zipfile containing cAmaretto results
}
\description{
cAMARETTO_ExportResults
}
\examples{
try(
cAMARETTO_ExportResults(  cAMARETTOresults,  cAMARETTOnetworkM,
  cAMARETTOnetworkC,   output_address="./")
)
}
