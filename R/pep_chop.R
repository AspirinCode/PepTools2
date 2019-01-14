#' Chop a vector of peptides into k-mers
#'
#' @param x A vector of peptides
#' @param k length of k-mers
#' @return A list of the possible k-mers generated for each peptide
#'
#' @examples
#'
#' peps = pep_ran(k = 9, n = 10)
#' kmers = pep_chop(x = peps, k = 3)
#'
#' @export
pep_chop = function(x, k){
  out = lapply(X = 1:length(x), FUN = function(i){
    n = nchar(x = x[i]) - k + 1
    current = sapply(X = 1:n, FUN = function(j){
      return( substr(x = x[i], start = j, stop = j + k - 1) )
    })
    return( current )
  })
  return(out)
}
