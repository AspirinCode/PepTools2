#' Mutate residues along peptide positions
#'
#' This is equivalent to making an alanine scan, where a given peptide is
#' substituted positionally with 'A'. Note, unique peptides are returned.
#'
#' @param pep A string representing a peptide
#'
#' @param pos Which position(s) to mutate
#'
#' @param res Which amino acid(s) to mutate to
#'
#' @param keep_wildtype Boolean should \code{pep} be included in the output
#'
#' @return A character vector with number of elements corresponding to the
#' number of peptides required to cover all mutations
#'
#' @examples
#' pep_mutate(pep = "RQGQDHPTM", pos = seq(1,9), res = AMINOACIDS$one, keep_wt = FALSE)
#' length(pep_mutate(pep = "RQGQDHPTM", pos = seq(1,9), res = AMINOACIDS$one, keep_wt = FALSE))
#' pep_mutate(pep = "RQGQDHPTM", pos = seq(1,9), res = AMINOACIDS$one, keep_wt = TRUE)
#' length(pep_mutate(pep = "RQGQDHPTM", pos = seq(1,9), res = AMINOACIDS$one, keep_wt = TRUE))
#'
#' @export
pep_mutate = function(pep, pos, res, keep_wt = FALSE){

  # Check arguments
  if( length(pep) > 1 ){
    stop("Only mutate one peptide at a time")
  }
  if( max(pos) > nchar(pep) ){
    stop("Trying to mutate positions beyond peptide length")
  }
  pep_check(pep = pep)

  # Do mutation
  p_mat = pep_split(pep = pep)
  res = unlist(strsplit(x = paste(res, collapse=''), split = ''))
  mutants = sapply(1:nrow(p_mat),function(i){
    sapply(1:length(pos), function(j){
      sapply(1:length(res), function(k){
        mut_ijk = p_mat[i,]
        mut_ijk[pos[j]] = res[k]
        mut_ijk = paste(mut_ijk, collapse = '')
        return(mut_ijk)
      })
    })
  })
  mutants = as.vector(mutants)

  # Add or remove wildtype peptide
  if( keep_wt ){
    mutants = c(pep, mutants)
  } else {
    mutants = mutants[mutants != pep]
  }

  # Check if output is empty
  if( length(mutants) == 0 ){
    stop("No peptides could be returned, e.g. pep_mutate('AAA', 1:3, 'A', keep_wt=FALSE)")
  }

  # Avoid duplicates
  mutants = unique(mutants)

  # Done, return!
  return(mutants)
}
