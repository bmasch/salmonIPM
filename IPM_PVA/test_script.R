# Load data
fish_data <- read.table(file.path("~", "SalmonIPM", "IPM_PVA", "fish_data.txt"), sep = "\t", header = T)

# Impute one NA value of S_tot_obs in Chamberlain 1986
fish_data$S_tot_obs[fish_data$pop == "Chamberlain" & fish_data$year == 1986] <- 
  mean(fish_data$S_tot_obs[fish_data$pop == "Chamberlain"][1:5], na.rm = T)

# Fit hierarchical spawner-recruit model to run reconstruction
RR_fit <- salmonIPM(fish_data = fish_data, model = "RR", chains = 3, iter = 1000, warmup = 500)

print(test_fit, pars = c("phi","R_hat"), include = F)

# Fit hierarchical IPM
IPM_fit <- salmonIPM(fish_data = fish_data, model = "IPM", chains = 3, iter = 2000, warmup = 1000, thin = 1, 
                     cores = 3, control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

print(IPM_fit, pars = c("mu_log_a","sigma_log_a",
                        "mu_log_b","sigma_log_b",
                        "sigma_log_phi","rho_log_phi",
                        "mu_p","sigma_alr_p","gamma_alr_p",
                        "mu_tau_alr_p","sigma_log_tau_alr_p","tau_alr_p",
                        "mu_sigma_proc","sigma_log_sigma_proc"))
launch_shinystan(IPM_fit)

# Look at R_tot process error in relation to spawner age composition




# Look at S_H_tot in relation to lagged brood year productivity anomalies



# Distribution of R_max and R_max*A