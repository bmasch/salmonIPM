#' Generate initial values for fitting either integrated or run-reconstruction spawner-recruit models in Stan.
#'
#' @param data Named list of input data for fitting either an integrated or run-reconstruction spawner-recruit model in \code{stan}, as returned by \code{stan_data}.
#' @param chains A positive integer specifying the number of Markov chains.
#' @param model One of \code{"IPM"} or \code{"RR"}, indicating whether to fit an integrated or run-reconstruction model.
#' @param pool_pops Logical, with default \code{TRUE}, indicating whether or not to treat the different populations as hierarchical rather than fixed/independent.
#' 
#' @importFrom stats aggregate na.omit
#' 
#' @return A list with initial starting values for all of the parameters and states in the Stan model.
#' 
#' @export
stan_init <- function(data, chains, model, pool_pops = TRUE) 
{
  if(model == "IPM")
  {
    with(data, {
      # for(i in names(fish_data)) assign(i, fish_data[[i]])
      S_tot_obs_noNA <- S_tot_obs
      if(N_S_obs < N)
        S_tot_obs[-which_S_obs] <- NA
      p_HOS_obs <- pmin(pmax(n_H_obs/(n_H_obs + n_W_obs), 0.001), 0.999)
      p_HOS_all <- rep(0,N)
      if(N_H > 0)
        p_HOS_all[which_H] <- p_HOS_obs
      min_age <- max_age - N_age + 1
      ages <- min_age:max_age
      q_obs <- sweep(n_age_obs, 1, rowSums(n_age_obs), "/")
      q_obs_NA <- apply(is.na(q_obs), 1, any)
      q_obs[q_obs_NA,] <- rep(colMeans(na.omit(q_obs)), each = sum(q_obs_NA))
      R <- matrix(NA, N, N_age)
      S_W_tot_obs <- S_tot_obs*(1 - p_HOS_all)
      B_rate_all <- rep(0,N)
      B_rate <- B_take_obs/(S_W_tot_obs[which_B]*(1 - q_obs[which_B,1]) + B_take_obs)
      B_rate[is.na(B_rate)] <- 0.1
      B_rate_all[which_B] <- B_rate
      
      # Maybe figure out a way to do this with a call to run_recon?
      for(i in 1:N)
        for(j in 1:N_age)
        {
          if(year[i] + ages[j] <= max(year[pop==pop[i]]))
          {
            b <- ifelse(j==1, 0, B_rate_all[i+ages[j]])
            f <- ifelse(j==1, 0, F_rate[i+ages[j]])
            R[i,j] <- S_W_tot_obs[i+ages[j]]*q_obs[i+ages[j],j]/((1 - b)*(1 - f))
          }
        }
      
      R <- pmax(R, min(1, R[R > 0], na.rm = T))
      R_tot <- rowSums(R)
      R_tot[is.na(R_tot)] <- max(R_tot, na.rm = T)
      p <- sweep(R, 1, R_tot, "/")
      p_NA <- apply(is.na(p), 1, any)
      p[p_NA, ] <- rep(colMeans(na.omit(p)), each = sum(p_NA))
      alr_p <- sweep(log(p[,1:(N_age-1)]), 1, log(p[,N_age]), "-")
      alr_p_z <- apply(alr_p, 2, scale)
      gamma_p <- aggregate(alr_p, list(pop), mean)
      gamma_alr_p <- aggregate(alr_p, list(pop), mean)[,-1]
      gamma_alr_p_z <- apply(gamma_alr_p, 2, scale)
      
      if(pool_pops)
      {
        return(lapply(1:chains, function(i)
          list(mu_log_a = runif(1,1,3),
               sigma_log_a = runif(1,0.1,0.5),
               log_a_z = array(runif(max(pop),-1,1), dim = max(pop)),
               mu_log_b = runif(1,-1,0),
               sigma_log_b = runif(1,0.1,0.5),
               log_b_z = array(runif(max(pop),-1,1), dim = max(pop)),
               beta_log_phi = array(rnorm(N_X,0,1), dim = N_X),
               rho_log_phi = runif(1,0.1,0.7),
               sigma_log_phi = runif(1,0.1,0.5),
               log_phi_z = array(rnorm(max(year),0,0.1), dim = max(year)),
               mu_sigma_proc = runif(1,0.05,2), sigma_log_sigma_proc = runif(1,0.1,0.5),
               log_sigma_proc_z = array(runif(max(pop),-1,1), dim = max(pop)),
               sigma_obs = runif(1,0.5,1),
               mu_p = colMeans(p), sigma_alr_p = array(runif(N_age-1,0.5,1), dim = N_age-1),
               gamma_alr_p_z = gamma_alr_p_z,
               mu_tau_alr_p = array(runif(N_age-1,0.5,1), dim = N_age-1),
               sigma_log_tau_alr_p = array(runif(N_age-1,0.5,1), dim = N_age-1),
               log_tau_alr_p_z = matrix(rnorm(max(pop)*(N_age-1),0,1), nrow = max(pop)),
               alr_p_z = alr_p_z,
               S_tot_init = rep(max(S_tot_obs_noNA), max_age*max(pop)),
               q_init = matrix(colMeans(q_obs), max_age*max(pop), 3, byrow = T),
               p_HOS = p_HOS_obs,
               log_R_tot_z = as.vector(scale(log(R_tot)))*0.1,
               B_rate = B_rate)))
      } else {
        return(lapply(1:chains, function(i)
          list(a = array(exp(runif(max(pop),1,3)), dim = max(pop)),
               b = array(exp(runif(max(pop),-1,0)), dim = max(pop)),
               beta_proc = matrix(rnorm(N_X*max(pop),0,1), max(pop), N_X),
               # logit_rho_proc = array(runif(max(pop),0,2), dim = max(pop)),
               rho_proc = array(runif(max(pop),0.1,0.7), dim = max(pop)),
               sigma_proc = array(runif(max(pop),0.05,2), dim = max(pop)), 
               sigma_obs = array(runif(max(pop),0.5,1), dim = max(pop)),
               gamma_p_arr = gamma_p,
               tau_alr_p = matrix(runif(max(pop)*(N_age-1),0.5,1), max(pop), N_age-1),
               alr_p_z = alr_p_z,
               S_tot_init = rep(max(S_tot_obs_noNA), max_age*max(pop)),
               q_init = matrix(colMeans(q_obs), max_age*max(pop), 3, byrow = T),
               p_HOS = p_HOS_obs,
               log_R_tot_z = as.vector(scale(log(R_tot)))*0.1,
               B_rate = B_rate)))
      }
    })
  } else if(model == "RR") {
    with(data, {
      if(pool_pops)
      {
        # This is currently not based on the input data
        return(lapply(1:chains, function(i)
          list(mu_log_a = runif(1,3,6), 
               sigma_log_a = runif(1,0.1,0.5),
               log_a_z = array(runif(max(pop),-1,1), dim = max(pop)), 
               mu_log_b = runif(1,-1,0),
               sigma_log_b = runif(1,0.1,0.5),
               log_b_z = array(runif(max(pop),-1,1), dim = max(pop)), 
               rho_log_phi = runif(1,0.1,0.7),
               sigma_log_phi = runif(1,0.1,0.5), 
               log_phi_z = array(rnorm(max(year),0,0.1), dim = max(year)),
               sigma = runif(1,0.1,2))))
      } else {
        return(lapply(1:chains, function(i)
          list(a = array(exp(runif(max(pop),1,3)), dim = max(pop)),
               b = array(exp(runif(max(pop),-1,0)), dim = max(pop)),
               rho = array(runif(max(pop),0.1,0.7), dim = max(pop)),
               sigma = array(runif(max(pop),0.5,1), dim = max(pop)))))
      }
    })
  }
}
