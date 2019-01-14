#' Compute the Kullback-Leibler Divergence
#'
#' @param peps A character vector of equal length peptides
#'
#' @param aa_background_distribution One of 'Viruses', 'Archaea', 'Bacteria',
#' 'Eukaryota', 'All', where 'ALL' is default. See ?BGFREQS for more.
#'
#' @return A list with 3 elements: i. \code{KLD_pos}: Named numeric vector with
#' the computed positional Kullback-Leibler Divergence, ii. \code{KLD_mat}: A
#' numeric matrix with computed Kullback-Leibler Divergence contribution from
#' each residue at each position, iii. \code{KLD_logo_aa_height}: The computed
#' Kullback-Leibler Divergence for logo visualisation, such that
#' \code{X[i,j] = frequency_matrix[i,j] * KLD_pos[i] * sign(KLD_mat[i,j]}
#'
#' @examples
#' kullback_leibler(peps = PEPTIDES, aa_background_distribution = 'All')
#'
#' @export
kullback_leibler = function(peps, aa_background_distribution = 'All'){

  # Check function arguments
  if( !(any(aa_background_distribution %in% rownames(BGFREQS))) ){
    stop("aa_background_distribution must be one of 'Viruses', 'Archaea',
         'Bacteria', 'Eukaryota', 'All'")
  }
  pep_check(pep = peps)

  # Compute frequency matrix
  frequency_matrix = pssm_freqs(pep = peps)

  # Set amino acid background distribution
  aa_bg_distn = unlist(BGFREQS[aa_background_distribution,])

  # Calculate PSSM of KLD values per residue per position
  KLD_mat = frequency_matrix
  KLD_mat[is.numeric(KLD_mat)] = NA
  for( i in 1:nrow(KLD_mat) ){
    for( j in 1:ncol(KLD_mat) ){
      res = colnames(KLD_mat)[j]
      p = frequency_matrix[i,j]
      q = aa_bg_distn[res]
      if( p == 0 ){
        KLD_mat[i,j] = 0
      } else {
        KLD_mat[i,j] = p * log2( p/q )
      }
    }
  }

  # Calculate positional KLD
  KLD_pos = rowSums(KLD_mat)

  # Calculate heights of letters in logo plot
  KLD_logo_aa_height = frequency_matrix
  KLD_logo_aa_height[is.numeric(KLD_logo_aa_height)] = NA
  for( i in 1:nrow(KLD_logo_aa_height) ){
    for( j in 1:ncol(KLD_logo_aa_height) ){
      KLD_logo_aa_height[i,j] = frequency_matrix[i,j] * KLD_pos[i] * sign(KLD_mat[i,j])
    }
  }
  return(list(KLD_pos = KLD_pos, KLD_mat = KLD_mat, KLD_logo_aa_height = KLD_logo_aa_height))
}
