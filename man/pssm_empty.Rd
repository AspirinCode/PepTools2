% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pssm_empty.R
\name{pssm_empty}
\alias{pssm_empty}
\title{Create empty Position Specific Scoring Matrix}
\usage{
pssm_empty(k, long = FALSE)
}
\arguments{
\item{k}{The number of peptide positions}

\item{long}{Add support for unknown 'X' and insertion '-'}
}
\value{
A numeric matrix with NA entries with dimensions \code{k} rows
and 20 columns
}
\description{
Create an empty Position Specific Scoring Matrix with \code{k} number of
rows and 20 columns sorted according to \code{ARNDCQEGHILKMFPSTWYV}
}
\examples{
pssm_empty(k = 9)

}
