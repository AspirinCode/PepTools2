#' Compute Position Specific Counts Matrix
#'
#' From a vector of equal length peptides, calculate the counts matrix,
#' counting the absolute occurrence of all amino acids
#'
#' @param peps A character vector of equal length peptides
#' @return A numeric matrix with \code{length(pep)} rows and 20 columns
#' @examples
#'
#'
#' @export
pssm_counts = function(pep){

  # Check input
  pep_check(pep = pep, require_length = TRUE)

  # Convert input peptides to peptide matrix
  pep_mat = pep_split(pep)

  # Set empty PSSM
  count_mat = pssm_empty(k = nchar(pep[1]))

  # Fill in counts
  for( res in colnames(count_mat) ){
    count_mat[,res] = colSums(pep_mat == res)
  }

  # Done, return!
  return(count_mat)
}
