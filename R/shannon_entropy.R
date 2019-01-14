#' Compute the Shannon Entropy
#'
#' @param peps A character vector of equal length peptides
#'
#' @return A list with 3 elements: i. \code{H_pos}: Named numeric vector with
#' the computed positional Shannon Entropy, ii. \code{H_mat}: A
#' numeric matrix with computed Shannon Entropy contribution from
#' each residue at each position, iii. \code{H_logo_aa_height}: The computed
#' Shannon Entropy for logo visualisation, such that
#' \code{X[i,j] = frequency_matrix[i,j] * H_pos[i] * sign(H_mat[i,j]}
#'
#' @examples
#' shannon_entropy(peps = PEPTIDES)
#'
#' @export
shannon_entropy = function(peps){

  # Check function arguments
  pep_check(pep = peps)

  # Compute frequency matrix
  frequency_matrix = pssm_freqs(pep = peps)

  # Calculate PSSM of H values per residue per position
  H_mat = frequency_matrix
  H_mat[is.numeric(H_mat)] = NA
  for( i in 1:nrow(H_mat) ){
    for( j in 1:ncol(H_mat) ){
      res = colnames(H_mat)[j]
      p = frequency_matrix[i,j]
      if( p == 0 ){
        H_mat[i,j] = 0
      } else {
        H_mat[i,j] = - p * log2( p )
      }
    }
  }

  # Calculate positional H
  H_pos = rowSums(H_mat)

  # Calculate heights of letters in logo plot
  H_logo_aa_height = frequency_matrix
  H_logo_aa_height[is.numeric(H_logo_aa_height)] = NA
  for( i in 1:nrow(H_logo_aa_height) ){
    for( j in 1:ncol(H_logo_aa_height) ){
      H_logo_aa_height[i,j] = frequency_matrix[i,j] * H_pos[i] * sign(H_mat[i,j])
    }
  }
  return(list(H_pos = H_pos, H_mat = H_mat, H_logo_aa_height = H_logo_aa_height))
}
