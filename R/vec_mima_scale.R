#' Min-max scale vector
#'
#' Min-max scale vector, such that all values v lie in [0;1]. If x is a numeric
#' vector, then min-max scale is: ( x - min(x) ) / ( max(x) - min(x) )
#' @param x A numeric vector
#' @param na_rm Boolean: Should NAs be removed
#' @examples
#' vec_mima_scale(rnorm(20))
#' @export
vec_mima_scale = function(x, na_rm = FALSE){
  if( !(is.vector(x) & is.numeric(x)) ){
    stop("Input must be a numeric matrix")
  }
  x_scaled = ( x - min(x,na.rm=na_rm) ) / ( max(x,na.rm=na_rm) - min(x,na.rm=na_rm) )
  return(x_scaled)
}
