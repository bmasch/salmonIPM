#' Simulate data under adult-to-adult integrated population model
#'
#' \code{IPM_adult_sim} simulates initial values for parameters and states in Stan.
#'
#' @param pars Model parameters to be used for simulations.
#' @param pop Population ID.
#' @param year Calendar year.
#' @param X Covariate(s).
#' @param N_age The number of adult age classes.
#' @param max_age Oldest adult age class.
#' @param A Area of spawning habitat.
#' @param F_rate Harvest rate.
#' @param B_rate Broodstock take rate.
#' @param SR_func Character code for the type of stocl-recruit model to fit. At present, only 'BH' for "Beverton-Holt' is allowed.
#' @param n_age_tot_obs Number of adults for age comp.
#' @param n_HW_tot_obs Number of adults for H vs W origin.
#' 
#' @return A list with initial starting values for all of the parameters and states in the Stan model.
#' 
#' @importFrom stats rbinom rlnorm rmultinom rnorm runif
#' 
#' @export
IPM_adult_sim <- function(pars, pop, year, X = NULL, N_age, max_age, A, 
                          F_rate, B_rate, SR_func = "BH", n_age_tot_obs, n_HW_tot_obs)
{
  # spawner-recruit functions
  BH <- function(a, Rmax, S, A) 
  {
    R <- a*S/(A + a*S/Rmax)
    return(R)
  }
  
  N <- length(pop)                        # number of observations 
  N_pop <- max(pop)                       # number of populations
  ages <- (max_age - N_age + 1):max_age   # adult ages
  if(is.null(X)) X <- matrix(0, nrow = max(year), ncol = 1)
  
  with(pars, {
    # parameters
    Sigma_log_aRmax <- diag(c(sigma_log_a, sigma_log_Rmax)^2)
    Sigma_log_aRmax[1,2] <- rho_log_aRmax*sigma_log_a*sigma_log_Rmax
    Sigma_log_aRmax[2,1] <- Sigma_log_aRmax[1,2]
    aRmax <- exp(mvrnorm(N_pop, c(mu_log_a, mu_log_Rmax), Sigma_log_aRmax))
    a <- aRmax[,1]
    Rmax <- aRmax[,2]
    log_phi <- rep(0, max(year))
    log_phi[1] <- rnorm(1, 0, sigma_log_phi/sqrt(1 - rho_log_phi^2))
    for(i in 2:length(log_phi))
      log_phi[i] <- rnorm(1, rho_log_phi*log_phi[i-1], sigma_log_phi)
    phi <- exp(log_phi + X %*% beta_log_phi)
    mu_alr_p <- log(mu_p[1:(N_age-1)]) - log(mu_p[N_age])
    Sigma_gamma <-  diag(sigma_gamma^2) * (L_gamma %*% t(L_gamma))
    gamma <- mvrnorm(N_pop, mu_alr_p, Sigma_gamma)
    Sigma_alr_p <- diag(sigma_alr_p^2) * (L_alr_p %*% t(L_alr_p))
    alr_p <- t(apply(gamma[pop,], 1, function(x) mvrnorm(1, x, Sigma_alr_p)))
    e_alr_p <- exp(cbind(alr_p, 0))
    p <- sweep(e_alr_p, 1, rowSums(e_alr_p), "/")
    
    # Simulate recruits and calculate total spawners
    # and spawner age distributions
    S_W <- matrix(NA, N, N_age)        # true wild spawners by age  
    S_W_tot <- vector("numeric",N)     # true total wild spawners
    S_H_tot <- vector("numeric",N)     # true total hatchery spawners
    S_tot <- vector("numeric",N)       # true total spawners
    R_tot_hat <- vector("numeric",N)   # expected recruits
    R_tot <- vector("numeric",N)       # true recruits
    B_take <- vector("numeric",N)      # adult broodstock removals
    
    for(i in 1:N)
    {
      if(year[i] - min(year[pop==pop[i]]) <= max_age)
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
      S_H_tot[i] <- S_W_tot[i]*p_HOS[i]/(1 - p_HOS[i])
      S_tot[i] <- S_W_tot[i] + S_H_tot[i]
      R_tot_hat[i] <- switch(SR_func,
                             BH = A[i]*BH(a[pop[i]], Rmax[pop[i]], S_tot[i], A[i]))
      R_tot[i] <- rlnorm(1, log(R_tot_hat[i]), sigma_proc)*phi[year[i]]
    }
    
    S_tot_obs <- rlnorm(N, log(S_tot), sigma_obs)           # obs total spawners
    q <- sweep(S_W, 1, S_W_tot, "/")                        # true spawner age distn 
    n_age_tot_obs <- pmax(round(pmin(n_age_tot_obs, S_tot_obs)), 1)  # cap age samples at pop size
    n_HW_tot_obs <- pmax(round(pmin(n_HW_tot_obs, S_tot_obs)), 1)    # cap H/W samples at pop size
    n_age_obs <- t(sapply(1:N, function(i) rmultinom(1, n_age_tot_obs[i], q[i,]))) # obs wild age frequencies
    dimnames(n_age_obs)[[2]] <- paste0("n_age", ages, "_obs")
    n_H_obs <- rbinom(N, n_HW_tot_obs, p_HOS)               # obs count of hatchery spawners
    n_W_obs <- n_HW_tot_obs - n_H_obs                       # obs count of wild spawners
    
    return(list(sim_dat = data.frame(pop = pop, A = A, year = year, fit_p_HOS = p_HOS > 0,
                                     S_tot_obs = S_tot_obs, n_age_obs, 
                                     n_H_obs = n_H_obs, n_W_obs = n_W_obs, 
                                     B_take_obs = B_take, F_rate = F_rate),
                pars_out = c(pars, list(S_W = S_W, a = a, Rmax = Rmax, phi = phi, 
                                        p_HOS = p_HOS, p = p, R_tot_hat = R_tot_hat, R_tot = R_tot))))
  })
}
