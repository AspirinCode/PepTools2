% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aa_translate.R
\name{aa_translate}
\alias{aa_translate}
\title{Translate between amino acid residue symbols}
\usage{
aa_translate(x, to)
}
\arguments{
\item{x}{A character vector of valid amino acid residue symbols either in
one-, three- or full format. Valid amino acid residue symbols are
\code{ARNDCQEGHILKMFPSTWYV}}

\item{to}{Which symbol to translate to, one of
\code{'one', 'three', 'full', 'chemistry'}}
}
\value{
A character vector with the translated symbols
}
\description{
\code{aa_translate()} tranlates back and forth between 'one', 'three', 'full'
and 'chemistry' representations of the 20 proteogenic amino acids
}
\examples{
aa_translate(AMINOACIDS$full, to = 'one')
aa_translate(AMINOACIDS$full, to = 'three')
aa_translate(AMINOACIDS$full, to = 'chemistry')

}
