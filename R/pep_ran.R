################################################################################
################################################################################
################################################################################
#' Generate \code{n} random \code{k}-mers (peptides)
#'
#' @param n The number of peptide(s) to generate
#' @param k The length of the peptide(s)
#' @return A character vector of n peptides
#' @examples
#' pep_ran(n = 1000, k = 9)
#'
#' @export
pep_ran = function(n, k){

  # A bit tricky to avoid looping over millions of peptides

  # Define the standard 20 amino acid residue characters
  res_chars = AMINOACIDS$one

  # First we generate one long vector with all the samples residues
  smpl_chars = sample(x = res_chars, size = n*k, replace = TRUE)

  # Then we collapse into one long string
  smpl_string = paste0(smpl_chars, collapse = '')

  # Now we generate indices corresponding to extracting every 9 characters
  to_index   = seq(k, k*n, by = k)
  from_index = to_index - (k-1)

  # and extract using sub string
  kmer_peptides = substring(text = smpl_string, first = from_index, last = to_index)

  # Done! Return!
  return(kmer_peptides)
}
################################################################################
################################################################################
################################################################################
