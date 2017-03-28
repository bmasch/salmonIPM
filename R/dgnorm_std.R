#' Standardized density function for generalized normal (Subbotin) distribution
#'
#' Standardized version of \code{dgnorm}.
#' 
#' @param x variate
#' @param mu mean of the distribution
#' @param sigma standard deviation of the distribution
#' @param shape shape of the distribution
#' 
#' @export
dgnorm_std <- function(x, mu = 0, sigma = 1, shape = 2)
{
  dgnorm(x, mu, sigma, shape)/dgnorm(mu, mu, sigma, shape)
}
