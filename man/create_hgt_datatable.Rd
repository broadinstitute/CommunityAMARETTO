% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cAMARETTO_HTMLreport.R
\name{create_hgt_datatable}
\alias{create_hgt_datatable}
\title{create_hgt_datatable}
\usage{
create_hgt_datatable(output_hgt, com_table = FALSE, ComNr = 1)
}
\arguments{
\item{output_hgt}{hyper geoetric test table}

\item{com_table}{TRUE if it is for community page, FALSE if index page.}

\item{ComNr}{community number}
}
\value{
DataTable
}
\description{
create_hgt_datatable
}
\examples{
try(create_hgt_datatable(output_hgt, com_table=FALSE, ComNr = 1))
}
