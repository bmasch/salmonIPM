# Density functions for generalized normal (Subbotin) distribution
dgnorm <- function(x, mu = 0, sigma = 1, shape = 2)
{
  (shape/2*sigma*gamma(1/shape))*exp(-(abs(x - mu)/sigma)^shape)
}

dgnorm_std <- function(x, mu = 0, sigma = 1, shape = 2)
{
  dgnorm(x, mu, sigma, shape)/dgnorm(mu, mu, sigma, shape)
}
