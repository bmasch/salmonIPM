options(device = windows)
library(lme4)
library(rethinking)

#-----------------
# PLOTS
#-----------------

# Load data
fish_data <- read.table(file.path("~", "SalmonIPM", "IPM_PVA", "fish_data.txt"), sep = "\t", header = T)
fish_data <- fish_data[fish_data$year >= 1960,]

# Impute one NA value of S_tot_obs in Chamberlain 1986
fish_data$S_tot_obs[fish_data$pop == "Chamberlain" & fish_data$year == 1986] <- 
  mean(fish_data$S_tot_obs[fish_data$pop == "Chamberlain"][1:5], na.rm = T)

# Load table of intrinsic potential by pop,
# add weighted area to fish_data
IP <- read.table(file.path("~", "SalmonIPM", "IPM_PVA", "IP.txt"), sep = "\t", header = T)
fish_data$Aw <- IP$Aw[match(fish_data$code, IP$code)]

# Create alternate version where A = 1 for all pops
fish_data2 <- fish_data
fish_data2$A <- 1

# Load catch data by year
catch_data0 <- read.table(file.path("~", "SalmonIPM", "IPM_PVA", "catch_data.txt"), sep = "\t", header = T)
catch_data <- catch_data0[catch_data0$year %in% fish_data$year,]
catch_data <- catch_data[order(catch_data$year),]

#-----------------
# FIT MODELS
#-----------------

# Fit non-hierarchical spawner-recruit model to run reconstruction
RR_fix <- salmonIPM(fish_data = fish_data, model = "RR", pool_pops = FALSE, chains = 3, iter = 1000, warmup = 500,
                    control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))
print(RR_fix, pars = c("phi","R_hat"), include = FALSE)

windows(width = 14, height = 10)
par(mfrow = c(4,6))
for(i in 1:24)
{
  hist(extract1(RR_fix,"a")[,i], 20, prob = T, xlab = "a", ylab = "Probability density",
       main = levels(fish_data$pop)[i])
  curve(dnorm(x,0,20), add = T, col = "red")
}

windows(width = 14, height = 10)
par(mfrow = c(4,6))
for(i in 1:24)
  plot(extract1(RR_fix,"a")[,i], extract1(RR_fix,"Rmax")[,i], xlab = "a", ylab = "Rmax",
       main = levels(fish_data$pop)[i], log = "xy")

windows(width = 14, height = 10)
par(mfrow = c(4,6))
rr <- run_recon(fish_data)
which_fit <- stan_data(fish_data, model = "RR")$which_fit
R_hat <- matrix(NA, nrow(extract1(RR_fix,"R_hat")), nrow(fish_data))
R_hat[,which_fit] <- extract1(RR_fix,"R_hat")[,which_fit]
for(i in 1:24)
{
  plot(rr$S[as.numeric(rr$pop)==i], rr$R[as.numeric(rr$pop)==i],
       xlab = "S", ylab = "R", main = levels(fish_data$pop)[i], pch = "")
  cc <- "blue"
  cc <- col2rgb(cc)
  cc <- rgb(cc[1], cc[2], cc[3], alpha = 0.3*255, maxColorValue = 255)
  for(j in sample(100))
    lines(sort(rr$S[as.numeric(rr$pop)==i], na.last = T),
          R_hat[j,as.numeric(rr$pop)==i][order(rr$S[as.numeric(rr$pop)==i], na.last = T)], col = cc)
  points(rr$S[as.numeric(rr$pop)==i], rr$R[as.numeric(rr$pop)==i], pch = 16)
}
rm(rr);rm(R_hat);rm(which_fit);rm(cc)

# Fit hierarchical spawner-recruit model to run reconstruction
RR_fit <- salmonIPM(fish_data = fish_data, model = "RR", pool_pops = TRUE, chains = 3, iter = 1000, warmup = 500,
                    control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))
print(RR_fit, pars = c("phi","R_hat"), include = FALSE)

windows(width = 14, height = 10)
par(mfrow = c(4,6))
for(i in 1:24)
  plot(extract1(RR_fit,"a")[,i], extract1(RR_fit,"Rmax")[,i], xlab = "a", ylab = "Rmax",
       main = levels(fish_data$pop)[i], log = "xy")

windows(width = 14, height = 10)
par(mfrow = c(4,6))
rr <- run_recon(fish_data)
which_fit <- stan_data(fish_data, model = "RR")$which_fit
R_hat <- matrix(NA, nrow(extract1(RR_fit,"R_hat")), nrow(fish_data))
R_hat[,which_fit] <- extract1(RR_fit,"R_hat")[,which_fit]
for(i in 1:24)
{
  plot(rr$S[as.numeric(rr$pop)==i], rr$R[as.numeric(rr$pop)==i],
       xlab = "S", ylab = "R", main = levels(fish_data$pop)[i], pch = "")
  cc <- "blue"
  cc <- col2rgb(cc)
  cc <- rgb(cc[1], cc[2], cc[3], alpha = 0.3*255, maxColorValue = 255)
  for(j in sample(100))
    lines(sort(rr$S[as.numeric(rr$pop)==i], na.last = T),
          R_hat[j,as.numeric(rr$pop)==i][order(rr$S[as.numeric(rr$pop)==i], na.last = T)], col = cc)
  points(rr$S[as.numeric(rr$pop)==i], rr$R[as.numeric(rr$pop)==i], pch = 16)
}
rm(rr);rm(R_hat);rm(which_fit);rm(cc)

# Fit hierarchical IPM
IPM_fit <- salmonIPM(fish_data = fish_data, model = "IPM", 
                     pars = c("mu_log_a","sigma_log_a","a",
                              "mu_log_Rmax","sigma_log_Rmax","Rmax","rho_log_aRmax",
                              "sigma_log_phi","rho_log_phi","phi",
                              "mu_p","sigma_gamma","R_gamma","gamma",
                              "sigma_alr_p","R_alr_p","p","sigma_proc","sigma_obs",
                              "S_tot","S_W_tot","S_H_tot","R_tot_hat","R_tot"),
                     chains = 3, iter = 2000, warmup = 1000, 
                     control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

print(IPM_fit, pars = c("mu_log_a","sigma_log_a",
                        "mu_log_Rmax","sigma_log_Rmax","rho_log_aRmax",
                        "sigma_log_phi","rho_log_phi",
                        "mu_p","sigma_gamma","R_gamma","gamma",
                        "sigma_alr_p","R_alr_p","sigma_proc","sigma_obs"))
launch_shinystan(IPM_fit)

# Fit non-hierarchical IPM
IPM_fix <- salmonIPM(fish_data = fish_data, model = "IPM", pool_pops = FALSE,
                     pars = c("a","Rmax","sigma_proc","rho_proc","sigma_proc","sigma_obs",
                              "gamma","sigma_alr_p","R_alr_p","p",
                              "S_tot","R_tot_hat","R_tot","log_R_tot_proc"),
                     chains = 3, iter = 2000, warmup = 1000, 
                     control = list(adapt_delta = 0.95, max_treedepth = 12))

print(IPM_fix, pars = c("a","Rmax","sigma_proc","rho_proc","sigma_proc","sigma_obs",
                        "gamma","sigma_alr_p","R_alr_p"))
launch_shinystan(IPM_fix)

windows(width = 14, height = 10)
par(mfrow = c(4,6))
for(i in 1:24)
  plot(extract1(IPM_fix,paste0("a[",i,"]")), 
       extract1(IPM_fix,paste0("Rmax[",i,"]")),
       log = "xy", xlab = "a", ylab = "Rmax", main = levels(fish_data$pop)[i])

windows(width = 14, height = 10)
par(mfrow = c(4,6))
for(i in 1:24)
{
  mod <- IPM_fix
  A <- unique(fish_data$A[as.numeric(fish_data$pop)==i])
  plot(apply(extract1(mod,"R_tot")[,as.numeric(fish_data$pop)==i], 1, median), 
       A*extract1(mod,"Rmax")[,i],
       log = "xy", xlab = "median(R_tot)", ylab = "Rmax", main = levels(fish_data$pop)[i])
  abline(0,1,col="red")
  rm(A);rm(mod)
}

windows(width = 14, height = 10)
par(mfrow = c(4,6))
for(i in 1:24)
{
  mod <- IPM_fix
  A <- unique(fish_data$A[as.numeric(fish_data$pop)==i])
  S_tot <- extract1(mod,"S_tot")[,as.numeric(fish_data$pop)==i]/A
  R_tot <- extract1(mod,"R_tot")[,as.numeric(fish_data$pop)==i]/A
  R_tot_hat <- extract1(mod,"R_tot_hat")[,as.numeric(fish_data$pop)==i]/A
  plot(colMeans(S_tot), colMeans(R_tot),
       xlab = "S_tot/ha", ylab = "R_tot/ha", main = levels(fish_data$pop)[i], pch = "")
  cc <- "blue"
  cc <- col2rgb(cc)
  cc <- rgb(cc[1], cc[2], cc[3], alpha = 0.3*255, maxColorValue = 255)
  for(j in sample(100))
    lines(sort(S_tot[j,]), R_tot_hat[j,order(S_tot[j,])], col = cc)
  points(colMeans(S_tot), colMeans(R_tot), pch = 16)
  points(run_recon(fish_data)$S[as.numeric(fish_data$pop)==i]/A, 
         run_recon(fish_data)$R[as.numeric(fish_data$pop)==i]/A)
  rm(list=c("mod","S_tot","R_tot","R_tot_hat","A"))
}

windows()
pairs(sapply(c("a[13]","Rmax[13]","sigma_proc[13]","sigma_obs[13]"), function(v) extract1(IPM_fix,v)), log = "xy")


# Fit extended "harvest" IPM
IPM_F_fit <- salmonIPM(fish_data = fish_data, catch_data = catch_data, model = "IPM_F", 
                     pars = c("mu_log_a","sigma_log_a","a",
                              "mu_log_Rmax","sigma_log_Rmax","Rmax","rho_log_aRmax",
                              "sigma_log_phi","rho_log_phi","phi",
                              "mu_p","sigma_gamma","R_gamma","gamma",
                              "sigma_alr_p","R_alr_p","p","sigma_proc","sigma_obs",
                              "c1","c2","F_rate", "sigma_log_C",
                              "S_tot","S_W_tot","S_H_tot","R_tot_hat","R_tot"),
                     chains = 3, iter = 2000, warmup = 1000, 
                     control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

print(IPM_F_fit, pars = c("mu_log_a","sigma_log_a",
                        "mu_log_Rmax","sigma_log_Rmax","rho_log_aRmax",
                        "sigma_log_phi","rho_log_phi",
                        "mu_p","sigma_gamma","R_gamma","gamma","c1","c2","F_rate","sigma_log_C",
                        "sigma_alr_p","R_alr_p","sigma_proc","sigma_obs"))
launch_shinystan(IPM_F_fit)



#------------------------
# SIMULATE DATA AND FIT
#------------------------

### Use data simulator

dat <- stan_data(fish_data, model = "IPM")
N_sim <- 20  # number of simulated datasets
N_pop <- max(dat$pop)
# N_year <- 50
N <- dat$N
sim_pars <- list(mu_log_a = 2, sigma_log_a = 0.2,
                 mu_log_Rmax = 3, sigma_log_Rmax = 0.2,
                 rho_log_aRmax = 0.5,
                 beta_log_phi = 0, rho_log_phi = 0.6, sigma_log_phi = 0.2,
                 sigma_proc = 0.2, mu_p = c(0.1,0.5,0.4),
                 sigma_gamma = c(0.2,0.3), L_gamma = diag(2),
                 sigma_alr_p = c(0.2,0.2), L_alr_p = diag(2),
                 sigma_obs = 0.2,
                 p_HOS = replace(rep(0,dat$N), dat$which_H, 
                                 pmin(dat$n_H_obs/(dat$n_H_obs + dat$n_W_obs), 0.95)))

IPM_npp_fit_sim_pars <- list(median = NULL, CI.025 = NULL, CI.975 = NULL)
IPM_pp_fit_sim_pars <- list(median = NULL, CI.025 = NULL, CI.975 = NULL)
RR_npp_fit_sim_pars <- IPM_pp_fit_sim_pars
RR_pp_fit_sim_pars <- IPM_pp_fit_sim_pars
for(i in 1:3)
{
  RR_npp_fit_sim_pars[[i]] <- as.data.frame(cbind(matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("a[", 1:N_pop, "]"))),
                                                  matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("Rmax[", 1:N_pop, "]"))),
                                                  matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("rho[", 1:N_pop, "]"))), 
                                                  matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("sigma[", 1:N_pop, "]")))))
  RR_pp_fit_sim_pars[[i]] <- as.data.frame(cbind(mu_log_a = rep(NA,N_sim), sigma_log_a = NA,
                                                 matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("a[", 1:N_pop, "]"))),
                                                 mu_log_Rmax = NA, sigma_log_Rmax = NA, 
                                                 matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("Rmax[", 1:N_pop, "]"))),
                                                 rho_log_aRmax = NA, rho_log_phi = NA, sigma_log_phi = NA, sigma = NA))
  IPM_npp_fit_sim_pars[[i]] <- as.data.frame(matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("a[", 1:N_pop, "]"))),
                                             matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("Rmax[", 1:N_pop, "]"))),
                                             matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("rho_proc[", 1:N_pop, "]"))), 
                                             matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("sigma_proc[", 1:N_pop, "]"))), 
                                             matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("sigma_obs[", 1:N_pop, "]"))))
  IPM_pp_fit_sim_pars[[i]] <- as.data.frame(cbind(mu_log_a = rep(NA,N_sim), sigma_log_a = NA,
                                                  matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("a[", 1:N_pop, "]"))),
                                                  mu_log_Rmax = NA, sigma_log_Rmax = NA, 
                                                  matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("Rmax[", 1:N_pop, "]"))),
                                                  rho_log_aRmax = NA, rho_log_phi = NA, sigma_log_phi = NA, 
                                                  sigma_proc = NA, sigma_obs = NA))
}

for(i in 1:N_sim)
{
  cat("simulated dataset", i, "of", N_sim, "\n")
  
  sim_out <- IPM_adult_sim(pars = sim_pars,
                           pop = dat$pop,
                           year = dat$year,
                           N_age = dat$N_age,
                           max_age = dat$max_age,
                           A = dat$A,
                           F_rate = dat$F_rate,
                           B_rate = replace(rep(0,dat$N), dat$which_B, 
                                            dat$B_take_obs/(dat$S_tot_obs[dat$which_B] + dat$B_take_obs)),
                           n_age_tot_obs = rowSums(dat$n_age_obs),
                           n_HW_tot_obs = rep(0,N))
  sim_fish_data <- sim_out$sim_dat
  
  # Fit single-population spawner-recruit model and store estimates
  RR_npp_stan_sim <- salmonIPM(fish_data = sim_fish_data, model = "RR", pool_pops = FALSE,
                               chains = 3, iter = 1000, warmup = 500,
                               control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))
  
  RR_npp_fit_sim <- as.matrix(RR_npp_stan_sim, names(RR_npp_fit_sim_pars$median))
  RR_npp_fit_sim_pars$median[i,] <- apply(RR_npp_fit_sim, 2, median)
  RR_npp_fit_sim_pars$CI.025[i,] <- apply(RR_npp_fit_sim, 2, quantile, 0.025)
  RR_npp_fit_sim_pars$CI.975[i,] <- apply(RR_npp_fit_sim, 2, quantile, 0.975)
  
  # Fit hierarchical spawner-recruit model and store estimates
  RR_pp_stan_sim <- salmonIPM(fish_data = sim_fish_data, model = "RR", 
                              chains = 3, iter = 1000, warmup = 500,
                              control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))
  
  RR_pp_fit_sim <- as.matrix(RR_pp_stan_sim, names(RR_pp_fit_sim_pars$median))
  RR_pp_fit_sim_pars$median[i,] <- apply(RR_pp_fit_sim, 2, median)
  RR_pp_fit_sim_pars$CI.025[i,] <- apply(RR_pp_fit_sim, 2, quantile, 0.025)
  RR_pp_fit_sim_pars$CI.975[i,] <- apply(RR_pp_fit_sim, 2, quantile, 0.975)
  
  # Fit single-population IPM and store estimates
  IPM_npp_stan_sim <- salmonIPM(fish_data = sim_fish_data, model = "IPM", pool_pops = FALSE, 
                                chains = 3, iter = 1000, warmup = 500, 
                                control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))
  
  IPM_npp_fit_sim <- as.matrix(IPM_npp_stan_sim, names(IPM_npp_fit_sim_pars$median))
  IPM_npp_fit_sim_pars$median[i,] <- apply(IPM_npp_fit_sim, 2, median)
  IPM_npp_fit_sim_pars$CI.025[i,] <- apply(IPM_npp_fit_sim, 2, quantile, 0.025)
  IPM_npp_fit_sim_pars$CI.975[i,] <- apply(IPM_npp_fit_sim, 2, quantile, 0.975)
  
  # Fit hierarchical IPM and store estimates
  IPM_pp_stan_sim <- salmonIPM(fish_data = sim_fish_data, model = "IPM", 
                               chains = 3, iter = 1000, warmup = 500, 
                               control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))
  
  IPM_pp_fit_sim <- as.matrix(IPM_pp_stan_sim, names(IPM_pp_fit_sim_pars$median))
  IPM_pp_fit_sim_pars$median[i,] <- apply(IPM_pp_fit_sim, 2, median)
  IPM_pp_fit_sim_pars$CI.025[i,] <- apply(IPM_pp_fit_sim, 2, quantile, 0.025)
  IPM_pp_fit_sim_pars$CI.975[i,] <- apply(IPM_pp_fit_sim, 2, quantile, 0.975)
}

# Plot estimated and true values of hyperparameters for RR_pp vs IPM_pp
dev.new(width = 14, height = 7)
par(mfrow = c(2,4), mar = c(2,3,2,1), oma = c(3,2,0,0))
c1 <- "orangered3"
c1 <- col2rgb(c1)
c1 <- rgb(c1[1], c1[2], c1[3], maxColorValue = 255, alpha = 255*0.5)
c2 <- "blue4"
c2 <- col2rgb(c2)
c2 <- rgb(c2[1], c2[2], c2[3], maxColorValue = 255, alpha = 255*0.5)
for(i in c("mu_log_a","sigma_log_a","mu_log_Rmax","sigma_log_Rmax","rho_log_aRmax","rho_log_phi","sigma_log_phi"))
{
  theta <- sim_pars[[i]]
  RR_pp_median <- RR_pp_fit_sim_pars$median[,i]
  RR_pp_CI.025 <- RR_pp_fit_sim_pars$CI.025[,i]
  RR_pp_CI.975 <- RR_pp_fit_sim_pars$CI.975[,i]
  IPM_pp_median <- IPM_pp_fit_sim_pars$median[,i]
  IPM_pp_CI.025 <- IPM_pp_fit_sim_pars$CI.025[,i]
  IPM_pp_CI.975 <- IPM_pp_fit_sim_pars$CI.975[,i]
  
  plot(RR_pp_median, IPM_pp_median, pch = "", las = 1, cex.axis = 1.2, cex.lab = 1.5,
       xlim = range(RR_pp_CI.025, RR_pp_CI.975, IPM_pp_CI.025, IPM_pp_CI.975),
       ylim = range(RR_pp_CI.025, RR_pp_CI.975, IPM_pp_CI.025, IPM_pp_CI.975),
       xlab = "", ylab = "", main = i)
  abline(0, 1, col = "darkgray", lwd = 2)
  abline(v = theta, h = theta, col = "darkgray", lwd = 2)
  segments(RR_pp_CI.025, IPM_pp_median, x1 = RR_pp_CI.975, col = c1)
  segments(RR_pp_median, IPM_pp_CI.025, y1 = IPM_pp_CI.975, col = c2)
  points(RR_pp_median, IPM_pp_median, pch = 16, cex = 1.2)
  mtext("RR estimate", side = 1, line = 1.5, outer = T)
  mtext("IPM estimate", side = 2, outer = T)
}

rm(list = c("c1","c2","RR_pp_median","RR_pp_CI.025","RR_pp_CI.975",
            "IPM_pp_median","IPM_pp_CI.025","IPM_pp_CI.975","theta"))




### Use Stan model output to simulate new data and then fit

IPM_sim <- salmonIPM(fish_data = fish_data, model = "IPM", 
                     pars = c("mu_log_a","sigma_log_a","a",
                              "mu_log_Rmax","sigma_log_Rmax","Rmax","rho_log_aRmax",
                              "sigma_log_phi","rho_log_phi","phi",
                              "mu_p","sigma_gamma","R_gamma","gamma",
                              "sigma_alr_p","R_alr_p","p","q","sigma_proc","sigma_obs",
                              "S_tot","S_W_tot","S_H_tot","R_tot_hat","R_tot"),
                     chains = 3, iter = 200, warmup = 100, 
                     control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

N_sim <- 30  # number of simulated datasets
N_pop <- length(levels(fish_data$pop))
sim_pars <- as.data.frame(cbind(mu_log_a = rep(NA,N_sim), sigma_log_a = NA,
                                matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("a[", 1:N_pop, "]"))),
                                mu_log_Rmax = NA, sigma_log_Rmax = NA, 
                                matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("Rmax[", 1:N_pop, "]"))),
                                rho_log_aRmax = NA, rho_log_phi = NA, sigma_log_phi = NA, sigma_proc = NA, sigma_obs = NA,
                                matrix(NA, N_sim, 3, dimnames = list(NULL,c("mu_p[1]","mu_p[2]","mu_p[3]"))),
                                matrix(NA, N_sim, 2, dimnames = list(NULL,c("sigma_gamma[1]","sigma_gamma[2]"))),
                                matrix(NA, N_sim, 2, dimnames = list(NULL,c("sigma_alr_p[1]","sigma_alr_p[2]")))))

IPM_pp_fit_sim_pars <- list(median = NULL, CI.025 = NULL, CI.975 = NULL)
RR_pp_fit_sim_pars <- IPM_pp_fit_sim_pars
for(i in 1:3)
{
  IPM_pp_fit_sim_pars[[i]] <- sim_pars
  RR_pp_fit_sim_pars[[i]] <- as.data.frame(cbind(mu_log_a = rep(NA,N_sim), sigma_log_a = NA,
                                                 matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("a[", 1:N_pop, "]"))),
                                                 mu_log_Rmax = NA, sigma_log_Rmax = NA, 
                                                 matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("Rmax[", 1:N_pop, "]"))),
                                                 rho_log_aRmax = NA, rho_log_phi = NA, sigma_log_phi = NA, sigma = NA))
}

for(i in 1:N_sim)
{
  cat("simulated dataset", i, "of", N_sim, "\n")
  
  sim_fish_data <- fish_data
  nn <- sample(IPM_sim@sim$iter,1)
  for(j in names(sim_pars))
    sim_pars[i,j] <- extract1(IPM_sim, j)[nn]
  sim_fish_data$S_tot_obs <- rlnorm(nrow(sim_fish_data), log(extract1(IPM_sim,"S_tot")[nn,]), extract1(IPM_sim,"sigma_obs")[nn])
  n_HW_tot_obs <- sim_fish_data$n_H_obs + sim_fish_data$n_W_obs
  sim_fish_data$n_H_obs <- rbinom(nrow(sim_fish_data), n_HW_tot_obs, 
                                  extract1(IPM_sim,"S_H_tot")[nn,]/(extract1(IPM_sim,"S_W_tot")[nn,] + extract1(IPM_sim,"S_H_tot")[nn,]))
  sim_fish_data$n_W_obs <- n_HW_tot_obs - sim_fish_data$n_H_obs
  n_age_obs <- rowSums(sim_fish_data[,c("n_age3_obs","n_age4_obs","n_age5_obs")])
  q_true <- extract1(IPM_sim,"q")[nn,,]
  for(j in 1:nrow(sim_fish_data))
    sim_fish_data[j,c("n_age3_obs","n_age4_obs","n_age5_obs")] <- rmultinom(1, n_age_obs[j], q_true[j,])
  
  # Fit hierarchical spawner-recruit model and store estimates
  RR_pp_fit_sim <- salmonIPM(fish_data = sim_fish_data, model = "RR", 
                             chains = 3, iter = 1000, warmup = 500,
                             control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))
  
  RR_pp_fit_sim <- as.matrix(RR_pp_fit_sim, names(RR_pp_fit_sim_pars$median))
  RR_pp_fit_sim_pars$median[i,] <- apply(RR_pp_fit_sim, 2, median)
  RR_pp_fit_sim_pars$CI.025[i,] <- apply(RR_pp_fit_sim, 2, quantile, 0.025)
  RR_pp_fit_sim_pars$CI.975[i,] <- apply(RR_pp_fit_sim, 2, quantile, 0.975)
  
  # Fit hierarchical IPM and store estimates
  IPM_pp_fit_sim <- salmonIPM(fish_data = sim_fish_data, model = "IPM", 
                              chains = 3, iter = 1000, warmup = 500, 
                              control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))
  
  IPM_pp_fit_sim <- as.matrix(IPM_pp_fit_sim, names(IPM_pp_fit_sim_pars$median))
  IPM_pp_fit_sim_pars$median[i,] <- apply(IPM_pp_fit_sim, 2, median)
  IPM_pp_fit_sim_pars$CI.025[i,] <- apply(IPM_pp_fit_sim, 2, quantile, 0.025)
  IPM_pp_fit_sim_pars$CI.975[i,] <- apply(IPM_pp_fit_sim, 2, quantile, 0.975)
}

# Plot estimated and true values of hyperparameters
dev.new(width = 14, height = 7)
par(mfrow = c(2,4), mar = c(2,3,2,1), oma = c(3,2,0,0))
c1 <- "orangered3"
c1 <- col2rgb(c1)
c1 <- rgb(c1[1], c1[2], c1[3], maxColorValue = 255, alpha = 255*0.5)
c2 <- "blue4"
c2 <- col2rgb(c2)
c2 <- rgb(c2[1], c2[2], c2[3], maxColorValue = 255, alpha = 255*0.5)
for(i in c("mu_log_a","sigma_log_a","mu_log_Rmax","sigma_log_Rmax","rho_log_aRmax","rho_log_phi","sigma_log_phi"))
{
  theta <- sim_pars[,i]
  RR_pp_median <- RR_pp_fit_sim_pars$median[,i]
  RR_pp_CI.025 <- RR_pp_fit_sim_pars$CI.025[,i]
  RR_pp_CI.975 <- RR_pp_fit_sim_pars$CI.975[,i]
  IPM_pp_median <- IPM_pp_fit_sim_pars$median[,i]
  IPM_pp_CI.025 <- IPM_pp_fit_sim_pars$CI.025[,i]
  IPM_pp_CI.975 <- IPM_pp_fit_sim_pars$CI.975[,i]
  
  plot(RR_pp_median, IPM_pp_median, pch = "", las = 1, cex.axis = 1.2, cex.lab = 1.5,
       xlim = range(RR_pp_CI.025, RR_pp_CI.975, IPM_pp_CI.025, IPM_pp_CI.975),
       ylim = range(RR_pp_CI.025, RR_pp_CI.975, IPM_pp_CI.025, IPM_pp_CI.975),
       xlab = "", ylab = "", main = i)
  abline(0, 1, col = "darkgray", lwd = 2)
  abline(v = theta, h = theta, col = "darkgray", lwd = 2)
  segments(RR_pp_CI.025, IPM_pp_median, x1 = RR_pp_CI.975, col = c1)
  segments(RR_pp_median, IPM_pp_CI.025, y1 = IPM_pp_CI.975, col = c2)
  points(RR_pp_median, IPM_pp_median, pch = 16, cex = 1.2)
  mtext("RR estimate", side = 1, line = 1.5, outer = T)
  mtext("IPM estimate", side = 2, outer = T)
}

rm(list = c("c1","c2","RR_pp_median","RR_pp_CI.025","RR_pp_CI.975",
            "IPM_pp_median","IPM_pp_CI.025","IPM_pp_CI.975","theta"))

# ##
# windows()
# plot(IPM_pp_fit_sim_pars$median$mu_log_a, 1:30, pch = 16,
#      xlim = range(IPM_pp_fit_sim_pars$CI.025$mu_log_a, IPM_pp_fit_sim_pars$CI.975$mu_log_a),
#      xlab = "mu_log_a", ylab = "data set")
# segments(x0 = IPM_pp_fit_sim_pars$CI.025$mu_log_a, y0 = 1:30, x1 = IPM_pp_fit_sim_pars$CI.975$mu_log_a)
# abline(v = extract1(IPM_sim,"mu_log_a")[nn])
# abline(v = mean(log(extract1(IPM_sim,"a")[nn,])), col = "gray")
# ##
# windows()
# plot(extract1(IPM_sim,"a")[nn,], colMeans(IPM_pp_fit_sim_pars$median[,paste0("a[",1:24,"]")]))
# abline(0,1)
# ##
# windows()
# plot(IPM_pp_fit_sim_pars$median$mu_log_Rmax, 1:30, pch = 16,
#      xlim = range(IPM_pp_fit_sim_pars$CI.025$mu_log_Rmax, IPM_pp_fit_sim_pars$CI.975$mu_log_Rmax),
#      xlab = "mu_log_Rmax", ylab = "data set")
# segments(x0 = IPM_pp_fit_sim_pars$CI.025$mu_log_Rmax, y0 = 1:30, x1 = IPM_pp_fit_sim_pars$CI.975$mu_log_Rmax)
# abline(v = extract1(IPM_sim,"mu_log_Rmax")[nn])
# ##
# windows()
# plot(extract1(IPM_sim,"Rmax")[nn,], colMeans(IPM_pp_fit_sim_pars$median[,paste0("Rmax[",1:24,"]")]))
# abline(0,1)
# ##



#-----------------
# PLOTS
#-----------------

# Look at R_tot process error in relation to spawner age composition

log_R_tot_z <- extract1(IPM_fit, "log_R_tot_z")
sigma_proc <- extract1(IPM_fit, "sigma_proc")[,stan_data(fish_data, model = "IPM")$pop]
log_R_tot_err <- colMeans(log_R_tot_z*sigma_proc)
qq <- matrix(stan_mean(IPM_fit, "q"), ncol = 3, byrow = T)
mean_age <- qq %*% matrix(3:5, ncol = 1)

dev.new(width = 10, height = 7)
par(mfrow = c(4,6), mar = c(2,2,3,1), oma = c(3,3,0,0))
for(i in levels(fish_data$pop))
{
  plot(mean_age[fish_data$pop==i], log_R_tot_err[fish_data$pop==i], xlab = "", ylab = "", main = i,
       las = 1, cex.lab = 1.2)
  abline(coef(lm(log_R_tot_err[fish_data$pop==i] ~ mean_age[fish_data$pop==i])))
}
mtext("Mean spawner age", side = 1, outer = TRUE, line = 1, cex = 1.5)
mtext("Recruitment process error residual", side = 2, outer = TRUE, line = 1, cex = 1.5)

summary(lmer(log_R_tot_err ~ mean_age + (mean_age | pop), data = fish_data))


# Look at S_H_tot in relation to lagged brood year productivity anomalies
# There's a decent positive relationship, esp. at lag 4, in some pops (others are a mess).
# Including pop-specific recruitment process errors does not improve the relationship.
# This could provide a non-mechanistic "process model" that could be used to simluate
# future hatchery spawner abundance.

dat <- stan_data(fish_data, model = "IPM")
year <- dat$year
L <- 4    # assume L is dominant age of hatchery spawners
year_lag <- year - L   
year_lag[year_lag <= 0] <- NA
S_H_tot <- stan_mean(IPM_fit, "S_H_tot")
log_phi_lag <- colMeans(log(extract1(IPM_fit, "phi")))[year_lag]
# log_R_tot_z <- extract1(IPM_fit, "log_R_tot_z")
# sigma_proc <- extract1(IPM_fit, "sigma_proc")[,dat$pop]
# log_R_tot_err <- colMeans(log_R_tot_z*sigma_proc)
# log_R_tot_err_lag <- rep(NA, dat$N)
# for(i in 1:dat$N)
#   log_R_tot_err_lag[i] <- ifelse(is.na(year_lag[i]), NA,
#                                  log_R_tot_err[dat$pop==dat$pop[i] & dat$year==year_lag[i]])

dev.new(width = 10, height = 7)
par(mfrow = c(3,5), mar = c(2,2,3,1), oma = c(3,3,0,0))
for(i in dat$which_pop_H)
  plot(log_phi_lag[dat$pop==i & fish_data$fit_p_HOS], S_H_tot[dat$pop==i & fish_data$fit_p_HOS], 
       xlab = "", ylab = "", main = levels(fish_data$pop)[i], las = 1, cex.lab = 1.2, log = "y")

mtext(paste0("Brood year log productivity anomaly (lag ", L, ")"), side = 1, outer = TRUE, line = 1, cex = 1.5)
mtext("Hatchery spawner abundance", side = 2, outer = TRUE, line = 1, cex = 1.5)


# Distribution of b vs b/A and R_max vs R_max*A

dat <- stan_data(fish_data2, model = "IPM")
log_Rmax <- log(extract1(IPM_fit2,"Rmax"))

dev.new(width = 7, height = 7)
par(mfrow = c(2,1))

qqnorm(colMeans(log_Rmax), main = "log(Rmax [spawners/ha])")
qqline(colMeans(log_Rmax))
qqnorm(colMeans(log_Rmax)*tapply(dat$A, dat$pop, mean), main = "log(Rmax [spawners])")
qqline(colMeans(log_Rmax)*tapply(dat$A, dat$pop, mean))


