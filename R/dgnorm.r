#' Density function for generalized normal (Subbotin) distribution
#'
#' Density function for a generalized normal distribution.
#' 
#' @param x variate
#' @param mu mean of the distribution
#' @param sigma standard deviation of the distribution
#' @param shape shape of the distribution
#' 
#' @export
dgnorm <- function(x, mu = 0, sigma = 1, shape = 2)
{
  (shape/2*sigma*gamma(1/shape))*exp(-(abs(x - mu)/sigma)^shape)
}
