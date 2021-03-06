% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pep_ran.R
\name{pep_ran}
\alias{pep_ran}
\title{Generate \code{n} random \code{k}-mers (peptides)}
\usage{
pep_ran(n, k)
}
\arguments{
\item{n}{The number of peptide(s) to generate}

\item{k}{The length of the peptide(s)}
}
\value{
A character vector of n peptides
}
\description{
Generate \code{n} random \code{k}-mers (peptides)
}
\examples{
pep_ran(n = 1000, k = 9)

}
