#' Extend peptides in a vector of peptides
#'
#' Given a character vector of peptides of unequal length, pep_extend will add
#' 'X to the end (C-terminus) of the peptides, such that any peptide shorter
#' than the longest peptide will be extended to match the length
#'
#' @param x A character vector of peptides to be extended
#' @return A character vector of peptides of equal length
#' @examples
#' peps = c('NWAWK','TATIDISM','KGLFEPWRCRC','MNSMEEMRESG',
#'          'HPQRHCHAQWYEEG','KIHSELVSVSVWA','WAGDIHHVC','NHSTGPKEGGFIMY',
#'          'GGFDVHIRWHTYKCP','ITHSPNMLSMVDT')
#' nchar(peps)
#' peps_ext = pep_extend(peps)
#' nchar(peps_ext)
#'
#' @export
pep_extend = function(pep, k = NULL, ext_char = 'X'){
  pep_check(pep = pep, require_length = FALSE)
  if( !is.null(k) && k < max(nchar(pep)) ){
    stop("k must larger than or equal to the longest peptide")
  }
  if( nchar(ext_char) > 1 ){
    stop("ext_char must be a single character")
  }
  if( is.null(k) ){
    k = max(nchar(pep))
  }
  extended_peps = sapply(1:length(pep), function(i){
    k_extend = k - nchar(pep[i])
    paste0(pep[i], paste( rep(ext_char, k_extend), collapse = ''))
  })
  return(extended_peps)
}
