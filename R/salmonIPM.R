#' Fits an integrated or run-reconstruction spawner-recruit model.
#'
#' @param fish_data Data frame that includes the following \code{colnames}, in no particular order except where noted:
#' \describe{
#' \item{\code{pop}}{Numeric or character population ID.}
#' \item{\code{year}}{Numeric variable giving the year the fish spawned (i.e., the brood year).}
#' \item{\code{A}}{Spawning habitat size (either stream length or area). Will usually be time-invariant within a population, but need not be.}
#' \item{\code{S_tot_obs}}{Total number (not density) of wild and hatchery-origin spawners.}
#' \item{\code{n_age[min_age]_obs...n_age[max_age]_obs}}{Multiple columns of observed spawner age frequencies (i.e., counts), where [min_age] is the numeral age in years (total, not ocean age) of the youngest spawners.}
#' \item{\code{n_W_obs}}{Observed frequency of natural-origin spawners.}
#' \item{\code{n_H_obs}}{Observed frequency of hatchery-origin spawners.}
#' \item{\code{fit_p_HOS}}{Logical or 0/1 indicating for each row in fish_data whether the model should estimate p_HOS > 0. This is only required if model == "IPM".}
#' \item{\code{F_rate}}{Total harvest rate (proportion) of natural-origin fish.}
#' \item{\code{B_take_obs}}{Number of adults taken for hatchery broodstock.}
#' }
#' @param env_data Optional data frame whose variables are time-varying environmental covariates, sequentially ordered with each row corresponding to a unique year in fish_data.
#' @param catch_data Only if model == "IPM_F", a data frame with numeric columns
#' @param model One of \code{"IPM"}, \code{"RR"}, or \code{"IPM_F"}, indicating whether the data are intended for an integrated or run-reconstruction model or the integrated "harvest" model.
#' @param pool_pops Logical, with default \code{TRUE}, indicating whether or not to treat the different populations as hierarchical rather than fixed/independent. Must be TRUE if model == "IPM_F".
#' @param init A list of named lists of initial values to be passed to \code{stan}. If \code{NULL}, initial values will be automatically generated from the supplied data using \code{stan_init}. 
#' @param pars A vector of character strings specifying parameters to monitor. If NULL, default values are used. If a non-default value is supplied, it is the user's responsibility to make sure the parameters requested appear in the model configuration specified.
#' @param chains A positive integer specifying the number of Markov chains.
#' @param iter A positive integer specifying the number of iterations for each chain (including warmup).
#' @param warmup A positive integer specifying the number of warmup (aka burnin) iterations per chain. If step-size adaptation is on (which it is by default), this also controls the number of iterations for which adaptation is run (and hence these warmup samples should not be used for inference). The number of warmup iterations should not be larger than \code{iter}
#' @param thin A positive integer specifying the period for saving samples. The default is 1, which is usually the recommended value.
#' @param cores Number of cores to use when executing the chains in parallel. Defaults to 3.
#' @param ... Additional arguments to pass to \code{stan}. 
#' #' 
#' @return An object of class \code{stanfit} representing the fitted model. See \code{rstan::stan} for details.
#' 
#' @importFrom rstan stan
#'
#' @export
salmonIPM <- function(fish_data, fish_data_fwd = NULL, env_data = NULL, catch_data = NULL, model, pool_pops = TRUE, 
                      init = NULL, pars = NULL, chains, iter, warmup, thin = 1, cores = 3, ...)
{
  dat <- stan_data(fish_data, fish_data_fwd, env_data, catch_data, model)
  if(is.null(pars))
    pars <- switch(model, 
                   IPM = switch(ifelse(pool_pops, "Y", "N"),
                                Y = c("mu_log_a","sigma_log_a","a",
                                      "mu_log_Rmax","sigma_log_Rmax","Rmax","rho_log_aRmax",
                                      "beta_log_phi","sigma_log_phi","rho_log_phi","phi",
                                      "mu_p","sigma_gamma","R_gamma","gamma",
                                      "sigma_alr_p","R_alr_p","p",
                                      "p_HOS","B_rate_all",
                                      "sigma_proc","sigma_obs",
                                      "S_tot","R_tot","q"),
                                N = c("a","Rmax","beta_proc","rho_proc","sigma_proc",
                                      "gamma","sigma_alr_p","R_alr_p","p",
                                      "p_HOS","B_rate_all","sigma_obs","S_tot","R_tot","q")),
                   IPM_F = c("mu_log_a","sigma_log_a","a",
                                      "mu_log_Rmax","sigma_log_Rmax","Rmax","rho_log_aRmax",
                                      "beta_log_phi","sigma_log_phi","rho_log_phi","phi",
                                      "mu_p","sigma_gamma","R_gamma","gamma",
                                      "sigma_alr_p","R_alr_p","p",
                                      "p_HOS","c1","c2","F_rate","sigma_log_C","B_rate_all",
                                      "sigma_proc","sigma_obs",
                                      "S_tot","R_tot","q"),
                   RR = switch(ifelse(pool_pops, "Y", "N"),
                               Y = c("mu_log_a","sigma_log_a","a",
                                     "mu_log_Rmax","sigma_log_Rmax","Rmax","rho_log_aRmax",
                                     "rho_log_phi","sigma_log_phi","phi","sigma",
                                     "R_hat","S_sim","R_sim"),
                               N = c("a","Rmax","rho","sigma","R_hat","S_sim","R_sim")))
  
  stan_path <- file.path(path.package("salmonIPM"), "stan")
  
  fit <- stan(file = switch(model,
                            IPM = file.path(stan_path, ifelse(pool_pops, "IPM_adult_pp.stan", "IPM_adult_npp.stan")),
                            IPM_F = file.path(stan_path, ifelse(pool_pops, "IPM_adult_pp_F.stan", "IPM_adult_npp.stan")),
                            RR = file.path(stan_path, ifelse(pool_pops, "SR_RR_pp.stan", "SR_RR_npp.stan"))),
              data = dat, 
              init = stan_init(dat, chains, model, pool_pops), 
              pars = pars,
              chains = chains, iter = iter, warmup = warmup, thin = thin, 
              cores = cores, ...)
}