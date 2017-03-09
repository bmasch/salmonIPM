# Convenience functions to simplify extracting posterior means
# or single parameters from stanfit objects
stan_mean <- function(object, pars)
{
  mm <- get_posterior_mean(object, pars)
  return(mm[,ncol(mm)])
}
