#' Min-max scale matrix
#'
#' Min-max scale matrix, such that all values v lie in [0;1]. If X is a numeric
#' matrix, then min-max scale is: ( X - min(X) ) / ( max(X) - min(X) )
#' @param x A numeric matrix
#' @param na_rm Boolean: Should NAs be removed?
#' @examples
#' mat_mima_scale(matrix(rnorm(180),nrow=9,ncol=20))
#' @export
mat_mima_scale = function(x, na_rm = FALSE){
  if( !(is.matrix(x) & is.numeric(x)) ){
    stop("Input must be a numeric matrix")
  }
  x_scaled = ( x - min(x,na.rm=na_rm) ) / ( max(x,na.rm=na_rm) - min(x,na.rm=na_rm) )
  return(x_scaled)
}
