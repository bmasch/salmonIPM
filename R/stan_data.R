#' Assemble input data for integrated or run-reconstruction spawner-recruit models.
#'
#' @param fish_data Data frame that includes the following \code{colnames}, in no particular order except where noted:
#' \describe{
#' \item{\code{pop}}{Numeric or character population ID.}
#' \item{\code{year}}{Integer variable giving the year the fish spawned (i.e., the brood year).}
#' \item{\code{A}}{Spawning habitat size (either stream length or area). Will usually be time-invariant within a population, but need not be.}
#' \item{\code{S_tot_obs}}{Total number (not density) of wild and hatchery-origin spawners.}
#' \item{\code{n_age_minAge...n_age_maxAge}}{Multiple columns of observed spawner age frequencies (i.e., counts), where minAge (maxAge) is the numeral age in years (total, not ocean age) of the youngest (oldest) spawners.}
#' \item{\code{n_W_obs}}{Observed frequency of natural-origin spawners.}
#' \item{\code{n_H_obs}}{Observed frequency of hatchery-origin spawners.}
#' \item{\code{fit_p_HOS}}{Logical or 0/1 indicating for each row in fish_data whether the model should estimate p_HOS > 0. This is only required if model == "IPM".}
#' \item{\code{F_rate}}{Total harvest rate (proportion) of natural-origin fish, only if model != "IPM_F".}
#' \item{\code{B_take_obs}}{Number of adults taken for hatchery broodstock.}
#' }
#' @param env_data Optional data frame whose variables are time-varying environmental covariates, sequentially ordered with each row corresponding to a unique year in \code{fish_data} (and \{fish_data_fwd}, if not \{NULL}).
#' @param fish_data_fwd Only if model == "IPM", optional data frame with the following \code{colnames}, representing "forward" or "future" simulations:
#' \describe{
#' \item{\code{pop}}{Numeric or character population ID. All values must also appear in \code{fish_data$pop}.}
#' \item{\code{year}}{Integer variable giving the year the fish spawned (i.e., the brood year). For each population in \code{fish_data_fwd$pop}, the first year appearing in \code{fish_data_fwd$year} must be one greater than the last year appearing in \code{fish_data$year}, i.e., \code{min(fish_data_fwd$year[fish_data_fwd$pop==j]) == max(fish_data$year[fish_data$pop==j]) + 1}.}
#' \item{\code{A}}{Spawning habitat size (either stream length or area). Will usually be time-invariant within a population, but need not be.}
#' \item{\code{F_rate}}{Total harvest rate (proportion) of natural-origin fish.}
#' \item{\code{B_rate}}{Total broodstock removal rate (proportion) of natural-origin fish.}
#' \item{\code{p_HOS}}{Proportion of hatchery-origin spawners.}
#' }
#' @param catch_data Only if model == "IPM_F", a data frame with numeric columns
#' \describe{
#' \item{\code{year}}{Year for fishery data. Must be identical to \code{unique(fish_data$year)}.}
#' \item{\code{R_F_obs}}{Total recruits to the fishery.}
#' \item{\code{C_obs}}{Total catch.}
#' }
#' @param model One of \code{"IPM"}, \code{"RR"}, or \code{"IPM_F"}, indicating whether the data are intended for an integrated or run-reconstruction model or the integrated "harvest" model.
#' 
#' @return A named list that can be passed to \code{stan} as the \code{data} argument. 
#' 
#' @export

stan_data <- function(fish_data, fish_data_fwd = NULL, env_data = NULL, catch_data = NULL, model)
{
  fish_data <- as.data.frame(fish_data)
  
  if(!is.null(fish_data_fwd))
  {
    if(model != "IPM")
      stop("Argument fish_data_fwd can only be specified if model == 'IPM'.\n")
    N_fwd <- nrow(fish_data_fwd)
    fish_data_fwd <- as.data.frame(fish_data_fwd)
    fish_data_fwd$pop <- factor(fish_data_fwd$pop, levels = levels(factor(fish_data$pop)))
    if(any(!fish_data_fwd$pop %in% fish_data$pop))
      stop("All populations in fish_data_fwd must appear in fish_data.\n")
    year_check <- tapply(fish_data$year, fish_data$pop, max)
    year_check <- year_check[names(year_check) %in% fish_data_fwd$pop]
    year_fwd_check <- tapply(fish_data_fwd$year, fish_data_fwd$pop, max)
    year_fwd_check <- year_fwd_check[names(year_fwd_check) %in% fish_data_fwd$pop]
    if(any(year_fwd_check != year_check + 1))
      stop("First year in fish_data_fwd must equal 1 + last year in fish_data for each population.\n")
    fish_data_fwd$pop <- as.numeric(fish_data_fwd$pop)
    fish_data_fwd$year <- as.numeric(factor(fish_data_fwd$year, 
                                            levels = levels(factor(c(fish_data$year, fish_data_fwd$year)))))
  }
  
  if(is.null(fish_data_fwd))
  {
    N_fwd <- 0
    fish_data_fwd <- data.frame(pop = 1, year = 1, A = 0, F_rate = 0, B_rate = 0, p_HOS = 0)
  }
  
  fish_data$pop <- as.numeric(factor(fish_data$pop))
  fish_data$year <- as.numeric(factor(fish_data$year))
  fish_data$fit_p_HOS <- as.logical(fish_data$fit_p_HOS)
  
  if(is.null(env_data))
    env_data <- matrix(0, max(fish_data$year, fish_data_fwd$year))
  
  if(nrow(env_data) != max(fish_data$year, fish_data_fwd$year)) 
    stop("Length of environmental time series does not equal number of brood years.\n")
  
  if(any(is.na(env_data)))
    stop("Missing values not allowed in environmental covariates.\n")
  
  if(model == 'IPM_F') 
  {
    if(is.null(catch_data) | any(is.na(catch_data)))
      stop("Missing values not allowed in total run size and catch with model == 'IPM_F'.\n")
  }
  
  for(i in c("pop","year","A","fit_p_HOS","B_take_obs"))
    if(any(is.na(fish_data[,i])))
      stop(paste0("Missing values not allowed in fish_data$", i, "\n"))
  
    if(any(is.na(fish_data_fwd)))
      stop("Missing values not allowed in fish_data_fwd.\n")
  
  max_age <- max(as.numeric(substring(names(fish_data)[grep("n_age", names(fish_data))], 6, 6)))
  F_rate_check <- tapply(fish_data$F_rate, fish_data$pop, function(x) any(is.na(x[-c(1:max_age)])))
  if(any(F_rate_check))
    stop(paste0("Missing values not allowed in fish_data$F_rate except in years 1:max_age", i, "\n"))
  
  if(any(is.na(fish_data$n_W_obs) != is.na(fish_data$n_H_obs)))
    stop(paste("Conflicting NAs in n_W_obs and n_H_obs in rows", 
               which(is.na(fish_data$n_W_obs) != is.na(fish_data$n_H_obs))), "\n")
  
  age_NA_check <- is.na(fish_data[,grep("n_age", names(fish_data))])
  if(any(!rowSums(age_NA_check) %in% c(0, nrow(age_NA_check))))
    stop(paste("Conflicting NAs in age frequency data in rows", 
               which(!rowSums(age_NA_check) %in% c(0, nrow(age_NA_check))), "\n"))
  
  if(!model %in% c("IPM", "IPM_F", "RR"))
    stop("Model must be either 'IPM', 'IPM_F', or 'RR'. \n")
  
  if(model == "IPM")
  {
    with(fish_data, {  
      dat <- list(N = nrow(fish_data),
                  pop = pop, 
                  year = year,
                  N_X = ncol(env_data), 
                  X = as.matrix(env_data),
                  N_pop_H = length(unique(pop[fit_p_HOS])),
                  which_pop_H = array(unique(pop[fit_p_HOS]), dim = length(unique(pop[fit_p_HOS]))),
                  N_S_obs = sum(!is.na(S_tot_obs)),
                  which_S_obs = array(which(!is.na(S_tot_obs)), dim = sum(!is.na(S_tot_obs))),
                  S_tot_obs = replace(S_tot_obs, is.na(S_tot_obs) | S_tot_obs==0, 1),
                  N_age = sum(grepl("n_age", names(fish_data))), 
                  max_age = max_age,
                  n_age_obs = as.matrix(fish_data[,grep("n_age", names(fish_data))]),
                  N_H = sum(fit_p_HOS),
                  which_H = array(which(fit_p_HOS), dim = max(sum(fit_p_HOS), 1)),
                  n_W_obs = array(n_W_obs[fit_p_HOS], dim = max(sum(fit_p_HOS), 1)),
                  n_H_obs = array(n_H_obs[fit_p_HOS], dim = max(sum(fit_p_HOS), 1)),
                  A = A,
                  F_rate = replace(F_rate, is.na(F_rate), 0),
                  N_B = sum(B_take_obs > 0),
                  which_B = array(which(B_take_obs > 0), dim = max(sum(B_take_obs > 0), 1)),
                  B_take_obs = B_take_obs[B_take_obs > 0],
                  N_fwd = N_fwd,
                  pop_fwd = fish_data_fwd$pop,
                  year_fwd = fish_data_fwd$year,
                  A_fwd = fish_data_fwd$A,
                  B_rate_fwd = fish_data_fwd$B_rate,
                  F_rate_fwd = fish_data_fwd$F_rate,
                  p_HOS_fwd = fish_data_fwd$p_HOS)
      
      if(dat$N_pop_H == 0) dat$which_pop_H <- array(1, dim = 1)
      if(dat$N_H == 0) 
      {
        dat$which_H <- array(1, dim = 1)
        dat$n_W_obs <- array(1, dim = 1)
        dat$n_H_obs <- array(1, dim = 1)
      }
      if(dat$N_B == 0)
      {
        dat$which_B <- array(1, dim = 1)
        dat$B_take_obs <- array(0, dim = 1)
      }

      dat$n_W_obs[is.na(dat$n_W_obs)] <- 0
      dat$n_H_obs[is.na(dat$n_H_obs)] <- 0
      dat$n_age_obs[is.na(dat$n_age_obs)] <- 0
      
      return(dat)
    })
  } else if(model == "IPM_F")
  {
    with(fish_data, {  
      dat <- list(N = nrow(fish_data),
                  pop = pop, 
                  year = year,
                  N_X = ncol(env_data), 
                  X = as.matrix(env_data),
                  N_pop_H = length(unique(pop[fit_p_HOS])),
                  which_pop_H = array(unique(pop[fit_p_HOS]), dim = length(unique(pop[fit_p_HOS]))),
                  N_S_obs = sum(!is.na(S_tot_obs)),
                  which_S_obs = array(which(!is.na(S_tot_obs)), dim = sum(!is.na(S_tot_obs))),
                  S_tot_obs = replace(S_tot_obs, is.na(S_tot_obs) | S_tot_obs==0, 1),
                  N_age = sum(grepl("n_age", names(fish_data))), 
                  max_age = max(as.numeric(substring(names(fish_data)[grep("n_age", names(fish_data))], 6, 6))),
                  n_age_obs = as.matrix(fish_data[,grep("n_age", names(fish_data))]),
                  N_H = sum(fit_p_HOS),
                  which_H = array(which(fit_p_HOS), dim = max(sum(fit_p_HOS), 1)),
                  n_W_obs = array(n_W_obs[fit_p_HOS], dim = max(sum(fit_p_HOS), 1)),
                  n_H_obs = array(n_H_obs[fit_p_HOS], dim = max(sum(fit_p_HOS), 1)),
                  A = A,
                  R_F_obs = array(catch_data$R_F_obs, dim = nrow(catch_data)),
                  C_obs = array(catch_data$C_obs, dim = nrow(catch_data)),
                  N_B = sum(B_take_obs > 0),
                  which_B = array(which(B_take_obs > 0), dim = max(sum(B_take_obs > 0), 1)),
                  B_take_obs = B_take_obs[B_take_obs > 0])
      
      if(dat$N_pop_H == 0) dat$which_pop_H <- array(1, dim = 1)
      if(dat$N_H == 0) 
      {
        dat$which_H <- array(1, dim = 1)
        dat$n_W_obs <- array(1, dim = 1)
        dat$n_H_obs <- array(1, dim = 1)
      }
      if(dat$N_B == 0)
      {
        dat$which_B <- array(1, dim = 1)
        dat$B_take_obs <- array(0, dim = 1)
      }

      dat$n_W_obs[is.na(dat$n_W_obs)] <- 0
      dat$n_H_obs[is.na(dat$n_H_obs)] <- 0
      dat$n_age_obs[is.na(dat$n_age_obs)] <- 0
      
      return(dat)
    })
  } else {
    recon_dat <- run_recon(fish_data)
    which_fit <- which(!is.na(recon_dat$R) & !is.na(recon_dat$S))
    N_fit <- length(which_fit)
    p <- aggregate(recon_dat[,grep("p_age", names(recon_dat))], list(pop = recon_dat$pop), mean, na.rm = TRUE)[,-1]
    
    with(recon_dat, {
      dat <- list(N = length(S),
                  pop = pop, 
                  year = year,
                  N_fit = N_fit,
                  which_fit = array(which_fit, dim = N_fit),
                  S = replace(S, S == 0 | is.na(S), 1),
                  R = replace(R, R == 0 | is.na(R), 1),
                  A = A,
                  S_NA = array(as.integer(is.na(S)), dim = length(S)),
                  R_NA = array(as.integer(is.na(R)), dim = length(R)),
                  N_age = sum(grepl("p_age", names(recon_dat))),
                  max_age = max(as.numeric(substring(names(recon_dat)[grep("p_age", names(recon_dat))], 6, 6))),
                  p = as.matrix(p))
      
      return(dat)
    })
  }
}
