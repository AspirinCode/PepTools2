#' Compute the Information Content
#'
#' @param peps A character vector of equal length peptides
#'
#' @return A list with 3 elements: i. \code{IC_pos}: Named numeric vector with
#' the computed positional Information Content, ii. \code{IC_mat}: A
#' numeric matrix with computed Information Content contribution from
#' each residue at each position, iii. \code{IC_logo_aa_height}: The computed
#' Information Content for logo visualisation, such that
#' \code{X[i,j] = frequency_matrix[i,j] * IC_pos[i] * sign(IC_mat[i,j]}
#'
#' @examples
#' information_content(peps = PEPTIDES)
#'
#' @export
information_content = function(peps){

  # Check function arguments
  pep_check(pep = peps)

  # Compute frequency matrix
  frequency_matrix = pssm_freqs(pep = peps)

  # Calculate PSSM of IC values per residue per position
  IC_mat = frequency_matrix
  IC_mat[is.numeric(IC_mat)] = NA
  for( i in 1:nrow(IC_mat) ){
    for( j in 1:ncol(IC_mat) ){
      res = colnames(IC_mat)[j]
      p = frequency_matrix[i,j]
      q = 1/20
      if( p == 0 ){
        IC_mat[i,j] = 0
      } else {
        IC_mat[i,j] = p * log2( p/q )
      }
    }
  }

  # Calculate positional IC
  IC_pos = rowSums(IC_mat)

  # Calculate heights of letters in logo plot
  IC_logo_aa_height = frequency_matrix
  IC_logo_aa_height[is.numeric(IC_logo_aa_height)] = NA
  for( i in 1:nrow(IC_logo_aa_height) ){
    for( j in 1:ncol(IC_logo_aa_height) ){
      IC_logo_aa_height[i,j] = frequency_matrix[i,j] * IC_pos[i] * sign(IC_mat[i,j])
    }
  }
  return(list(IC_pos = IC_pos, IC_mat = IC_mat, IC_logo_aa_height = IC_logo_aa_height))
}
