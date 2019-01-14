#' Check if a vector of peptides is valid.
#'
#' According to criteria:
#'
#'  i.   Is a vector with one or more elements
#'
#'  ii.  Is a character vector
#'
#'  iii. All elements in vector have the same number of characters (depending on the require_length argument)
#'
#'  iv.  Elements only contain allowed amino acid residues \code{'ARNDCQEGHILKMFPSTWYVX-'}
#'
#' @param pep A vector of peptides to be checked
#' @return TRUE if the vector passed the check, otherwise stop() is called
#' @examples
#' pep_check(c("RQGQDHPTM","RGQKTTDNA","NILYEYWDY"))
#' pep_check(c("RQGQDHPTM","RGQKTTDNA","NILYEYWDYX"), require_length = FALSE)
#' pep_check(c("RQGQDHPTM","RGQKTTDNA","NILYEYWDB"))
#'
#' @export
pep_check = function(pep, require_length = TRUE){
  # Check if 'pep' is a character vector
  if( !( is.vector(pep) & is.character(pep) ) ){
    stop("'pep' has to be a character vector with one or more peptides!")
  }
  # Check if all sequences are of equal length
  if( require_length & !all( nchar(pep) == nchar(pep[1])) ){
    stop("All peptides must have equal number of positions")
  }
  # Check if 'pep' contain non-allowed amino acid residue characters
  if( grepl(pattern = "[^ARNDCQEGHILKMFPSTWYVX-]", x = paste0(pep, collapse='')) ){
    stop("Non standard amino acid residue character found.
         Only 'ARNDCQEGHILKMFPSTWYVX-' allowed")
  }
  return(TRUE)
}
