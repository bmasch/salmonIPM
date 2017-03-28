#' Simulate data under adult-to-adult integrated population model
#'
#' \code{IPM_adult_sim} simulates initial values for parameters and states in Stan.
#'
#' @param \code{pars}{Model parameters to b eused for simulations.}
#' @param \code{pop}{Population ID.}
#' @param \code{year}{ }
#' @param \code{X}{ }
#' @param \code{N_age}{The number of adult age classes.}
#' @param \code{max_age}{Oldest adult age class.}
#' @param \code{S_H_tot}{ }
#' @param \code{A}{ }
#' @param \code{F_rate}{ }
#' @param \code{B_rate}{ }
#' @param \code{SR_func}{Character code for the type of stocl-recruit model to fit. At present, only 'BH' for "Beverton-Holt' is allowed.}
#' @param \code{n_age_tot_obs}{ }
#' @param \code{n_HW_tot_obs}{ }
#' 
#' @return A list with initial starting values for all of the parameters and states in the Stan model.
#' 
#' @export
IPM_adult_sim <- function(pars, pop, year, X, N_age, max_age, S_H_tot, A, 
                          F_rate, B_rate, SR_func = "BH", n_age_tot_obs, n_HW_tot_obs)
{
  # spawner-recruit functions
  BH <- function(a, b, S, A) 
  {
    R <- a*S/(A + b*S)
    return(R)
  }
  
  BH_Rmax <- function(a, b, A)
  {
    return(A*a/b)
  }
  
  N <- length(pop)                        # number of observations 
  N_pop <- max(pop)                       # number of populations
  which_pop_H <- unique(pop[S_H_tot > 0]) # populations with hatchery input
  ages <- (max_age - N_age + 1):max_age   # adult ages
  
  # parameters
  mu_log_a <- pars$mu_log_a
  sigma_log_a <- pars$sigma_log_a
  a <- rlnorm(N_pop, mu_log_a, sigma_log_a)
  mu_log_b <- pars$mu_log_b
  sigma_log_b <- pars$sigma_log_b
  b <- rlnorm(N_pop, mu_log_b, sigma_log_b)
  Rmax <- BH_Rmax(a, b, A)
  beta_log_phi <- pars$beta_log_phi
  rho_log_phi <- pars$rho_log_phi
  sigma_log_phi <- pars$sigma_log_phi
  mu_log_phi <- X %*% beta_log_phi
  log_phi <- rep(0, max(year))
  log_phi[1] <- rnorm(1, 0, sigma_log_phi/sqrt(1 - rho_log_phi^2))
  for(i in 2:max(year))
    log_phi[i] <- rnorm(1, rho_log_phi*log_phi[i-1], sigma_log_phi)
  phi <- exp(log_phi + X %*% beta_log_phi)
  mu_sigma_proc <- pars$mu_sigma_proc
  sigma_log_sigma_proc <- pars$sigma_log_sigma_proc
  sigma_proc <- rlnorm(N_pop, log(mu_sigma_proc), sigma_log_sigma_proc)
  sigma_obs <- pars$sigma_obs
  mu_alr_p <- log(pars$mu_p[1:(N_age-1)]) - log(pars$mu_p[N_age])
  sigma_alr_p <- pars$sigma_alr_p
  gamma_alr_p <- t(sapply(1:N_pop, function(i) rnorm(N_age - 1, mu_alr_p, sigma_alr_p)))
  mu_tau_alr_p <- pars$mu_tau_alr_p
  sigma_log_tau_alr_p <- pars$sigma_log_tau_alr_p
  tau_alr_p <- t(sapply(1:N_pop, function(i) rlnorm(N_age - 1, log(mu_tau_alr_p), sigma_log_tau_alr_p)))
  alr_p <- t(sapply(1:N, function(i) rnorm(N_age - 1, gamma_alr_p[pop[i],], tau_alr_p[pop[i],])))
  p <- sweep(cbind(exp(alr_p), 1), 1, rowSums(exp(alr_p)) + 1, "/")

  # Simulate recruits and calculate total spawners
  # and spawner age distributions
  S_W <- matrix(NA, N, 3)            # true wild spawners by age  
  S_W_tot <- vector("numeric",N)     # true total wild spawners
  S_tot <- vector("numeric",N)       # true total spawners
  R_tot_hat <- vector("numeric",N)   # expected recruits
  R_tot <- vector("numeric",N)       # true recruits
  B_take <- vector("numeric",N)      # adult broodstock removals
  
  for(i in 1:N)
  {
    if(year[i] <= max_age)
    {
      S_W[i,] <- rlnorm(N_age, log(Rmax[pop[i]]/N_age), 0.1) # initialize years 1:max_age
    } else
    {
      for(j in 1:N_age)
        S_W[i,j] <- R_tot[i-ages[j]]*p[i-ages[j],j]
      S_W[i,-1] <- S_W[i,-1]*(1 - F_rate[i])     # catch (assumes no take of age 1)
      B_take[i] <- B_rate[i]*sum(S_W[i,-1])
      S_W[i,-1] <- S_W[i,-1]*(1 - B_rate[i])     # broodstock removal (assumes no take of age 1)
    }
    S_W_tot[i] <- sum(S_W[i,])
    S_tot[i] <- S_W_tot[i] + S_H_tot[i]
    R_tot_hat[i] <- switch(SR_func,
                           BH = A*BH(a[pop[i]], b[pop[i]], S_tot[i], A)*phi[year[i]])
    R_tot[i] <- rlnorm(1, log(R_tot_hat[i]), sigma_proc[pop[i]])
  }
  
  S_tot_obs <- rlnorm(N, log(S_tot), sigma_obs)           # obs total spawners
  pHOS <- S_H_tot/(S_W_tot + S_H_tot)                     # true pHOS
  q <- sweep(S_W, 1, S_W_tot, "/")                        # true spawner age distn 
  n_age_tot_obs <- round(pmin(n_age_tot_obs, S_tot_obs))  # cap age samples at pop size
  n_HW_tot_obs <- round(pmin(n_HW_tot_obs, S_tot_obs))    # cap H/W samples at pop size
  n_age_obs <- t(sapply(1:N, function(i) rmultinom(1, n_age_tot_obs[i], q[i,]))) # obs wild age frequencies
  n_H_obs <- rbinom(N, n_HW_tot_obs, pHOS)                # obs count of hatchery spawners
  n_W_obs <- n_HW_tot_obs - n_H_obs                       # obs count of wild spawners
  
  return(list(sim_dat = list(N = N, pop = pop, year = year, X = X, S_tot_obs = S_tot_obs, 
                             n_age_obs = n_age_obs, n_H_obs = n_H_obs, n_W_obs = n_W_obs, 
                             F_rate = F_rate, B_take = B_take, A = A),
              pars_out = c(pars, list(S_W = S_W, a = a, b = b, phi = phi, 
                                      pHOS = pHOS, p = p, sigma_proc = sigma_proc,
                           R_tot_hat = R_tot_hat, R_tot = R_tot))))
}
