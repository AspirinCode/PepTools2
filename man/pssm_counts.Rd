% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pssm_counts.R
\name{pssm_counts}
\alias{pssm_counts}
\title{Compute Position Specific Counts Matrix}
\usage{
pssm_counts(pep)
}
\arguments{
\item{peps}{A character vector of equal length peptides}
}
\value{
A numeric matrix with \code{length(pep)} rows and 20 columns
}
\description{
From a vector of equal length peptides, calculate the counts matrix,
counting the absolute occurrence of all amino acids
}
\examples{


}
