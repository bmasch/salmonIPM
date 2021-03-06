#' Extract posterior means from Stan object
#'
#' Convenience function to simplify extracting posterior means or single
#' parameters from stanfit objects.
#' 
#' @param object A fitted Stan object.
#' @param pars Character vector with the names of parameters to summarize.
#' 
#' @return A scalar or vector of means of the posterior distribution of \code{pars}.
#' 
#' @importFrom rstan get_posterior_mean
#'
#' @export
stan_mean <- function(object, pars)
{
  mm <- get_posterior_mean(object, pars)
  return(mm[,ncol(mm)])
}
