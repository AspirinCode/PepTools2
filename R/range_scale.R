#' Linearly scales values to new range
#'
#' @param x A set of numbers as a numeric vector, matrix or data frame
#' @return The scaled values in same structure as was input
#' @examples
#' sample = rnorm(10)
#' range_scale(sample)
#'
#' @export
range_scale = function(x, new_min = -1, new_max = 1){
  b = (new_max - new_min) / (max(x) - min(x))
  a = new_max - b * max(x)
  y = a + b * x
  return(y)
}
