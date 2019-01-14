#' Compute frequency matrix
#'
#' From a vector of equal length peptides, calculate the corresponding frequency
#' matrix
#'
#' @param pep A character vector of equal length peptides
#'
#' @return A numeric matrix with \code{length(pep)} rows and 20 columns
#'
#' @examples
#' pssm_freqs(pep = pep_ran(n = 100, k = 9))
#'
#' @export
pssm_freqs = function(pep){

  # Check input arguments
  pep_check(pep)

  # Create count matrix
  count_mat = pssm_counts(pep = pep)

  # Compute frequency matrix
  f_mat = count_mat / rowSums(count_mat)

  # Done, return!
  return(f_mat)
}
