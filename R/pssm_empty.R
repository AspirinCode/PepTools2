#' Create empty Position Specific Scoring Matrix
#'
#' Create an empty Position Specific Scoring Matrix with \code{k} number of
#' rows and 20 columns sorted according to \code{ARNDCQEGHILKMFPSTWYV}
#'
#' @param k The number of peptide positions
#' @param long Add support for unknown 'X' and insertion '-'
#' @return A numeric matrix with NA entries with dimensions \code{k} rows
#' and 20 columns
#' @examples
#' pssm_empty(k = 9)
#'
#' @export
pssm_empty = function(k, long = FALSE){

  # Set aminoacids
  aa = AMINOACIDS$one

  # If long, add support for unknown 'X' and insertion '-'
  if( long ){ aa = c(aa, c('X', '-')) }
  return( matrix(data = NA, nrow = k, ncol = length(aa),
                 dimnames = list(paste0('p', 1:k), aa)) )
}

