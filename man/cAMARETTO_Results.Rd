% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cAMARETTO_Results.R
\name{cAMARETTO_Results}
\alias{cAMARETTO_Results}
\title{cAMARETTO_Results

To initiate the community AMARETTO this results functions performs
hyper geometric tests between all modules and gene sets.}
\usage{
cAMARETTO_Results(
  AMARETTOresults_all,
  NrCores = 1,
  output_dir = "./",
  gmt_filelist = NULL,
  drivers = FALSE
)
}
\arguments{
\item{AMARETTOresults_all}{A list of multiple
AMARETTO_Run outputs. The names are run names.}

\item{NrCores}{Nr of Cores that can be used to
calculated the results.}

\item{output_dir}{A directory that stores gmt files
for the runs}

\item{gmt_filelist}{A list with gmt files that are added into
the communities. The names of the list are used in the networks.
NULL if no list is added.}

\item{drivers}{Boolean that defines if only targets or drivers
and targets are used to calculate the HGT.}
}
\value{
a list with cAMARETTOresults data objects from multiple runs
}
\description{
cAMARETTO_Results

To initiate the community AMARETTO this results functions performs
hyper geometric tests between all modules and gene sets.
}
\examples{
try(
cAMARETTOresults <- cAMARETTO_Results(AMARETTOresults_all,
  gmt_filelist=list(ImmuneSignature = Cibersortgmt),
  NrCores = 4 ,
  output_dir = "./")
)
}
