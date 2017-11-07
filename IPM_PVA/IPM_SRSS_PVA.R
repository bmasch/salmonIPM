setwd(file.path("~","salmonIPM","IPM_PVA"))
options(device=windows)
library(salmonIPM)
library(corrplot)
library(magicaxis)
library(zoo)

#===========================================================================
# DATA
#===========================================================================

# Load data
fish_data <- read.table(file.path("~", "salmonIPM", "IPM_PVA", "fish_data.txt"), sep = "\t", header = T)
fish_data <- fish_data[!is.na(fish_data$B_take_obs),]
fish_data <- fish_data[order(fish_data$code, fish_data$year),]

# Load habitat area data and add area column (convert m2 to ha) to fish_data
IP <- read.table(file.path("~", "salmonIPM", "IPM_PVA", "IP.txt"), sep = "\t", header = T)
fish_data <- cbind(fish_data[,1:5], A = IP$A[match(fish_data$code, IP$code)]/1e4, fish_data[,-c(1:5)])

# Create dummy covariate for intervention analysis: pre/post-1970 (centered, not scaled)
pre_post_1970 <- data.frame(pre_post_1970 = as.numeric(sort(unique(fish_data$year)) >= 1970),
                       row.names = sort(unique(fish_data$year)))
pre_post_1970$pre_post_1970 <- pre_post_1970$pre_post_1970 - mean(pre_post_1970$pre_post_1970)

# Pad data with years through max_year
N_future_years <- 50
max_year <- max(fish_data$year) + N_future_years
year_aug <- sapply(tapply(fish_data$year, fish_data$pop, max), function(x) (x + 1):max_year)
pop_aug <- rep(names(year_aug), sapply(year_aug, length))
code_aug <- fish_data$code[match(pop_aug, fish_data$pop)]
MPG_aug <- fish_data$MPG[match(pop_aug, fish_data$pop)]
ESU_aug <- fish_data$ESU[match(pop_aug, fish_data$pop)]
A_aug <- rep(tapply(fish_data$A, fish_data$pop, mean), times = sapply(year_aug, length))
fish_data_aug <- data.frame(pop = pop_aug, code = code_aug, MPG = MPG_aug, ESU = ESU_aug, A = A_aug,
                            year = unlist(year_aug), type = "future", fit_p_HOS = 0, 
                            S_tot_obs = NA, n_age3_obs = 0, n_age4_obs = 0, n_age5_obs = 0,
                            n_W_obs = 0, n_H_obs = 0, p_HOS = 0, B_take_obs = 0, F_rate = 0,
                            row.names = NULL)
fish_data_aug <- rbind(cbind(type = "past", fish_data[,setdiff(names(fish_data_aug), "type")])[,names(fish_data_aug)], 
                       fish_data_aug)
fish_data_aug <- fish_data_aug[order(fish_data_aug$code, fish_data_aug$year),]
row.names(fish_data_aug) <- NULL

# Create "forward simulation" data
fish_data_fwd <- data.frame(pop = rep(pop_aug,2), year = rep(unlist(year_aug), 2),
                            A = rep(A_aug,2), F_rate = 0, B_rate = 0, p_HOS = 0)
row.names(fish_data_fwd) <- NULL


#--------------------------------------------------------
# Create data summary table 
#--------------------------------------------------------

Ph <- fish_data$p_HOS
S_tot_obs <- fish_data$S_tot_obs
pop <- fish_data$pop
table1 <- data.frame(Population=rep(NA,length(levels(pop))), Years=NA, Area=NA, St=NA, P.hatchery=NA)
table1$Population <- levels(pop)[order(fish_data$MPG[match(levels(pop), pop)])]
for(i in 1:nrow(table1))
{
  table1$Years[i] <- paste(min(fish_data$year[pop==table1$Population[i]]), "-", 
                           max(fish_data$year[pop==table1$Population[i]]), sep="")
  table1$St[i] <- paste0(round(median(S_tot_obs[pop==table1$Population[i]], na.rm = T),0), " (", 
                         paste0(round(quantile(S_tot_obs[pop==table1$Population[i]], c(0.05,0.95), na.rm = T),0), collapse = "-"),
                         ")")
  table1$P.hatchery[i] <- paste(round(mean(Ph[pop==table1$Population[i]]),2), 
                                ifelse(round(mean(Ph[pop==table1$Population[i]]),2)==0, "",
                                       paste0("(", 
                                              paste(round(quantile(Ph[pop==table1$Population[i]], c(0.05,0.95)),2), collapse = "-"),
                                              ")")))
}
table1$Area <- round(fish_data$A[match(table1$Population, pop)], 1)
rm(Ph);rm(S_tot_obs);rm(pop)
print(table1)
write.table(table1, "table1.txt", sep="\t", row.names=F)


#===========================================================================
# FIT RETROSPECTIVE MODELS
#
# Fit to observed data to check model behavior 
#===========================================================================

# Base model
set.seed(123)
IPM_pp1 <- salmonIPM(fish_data = fish_data, fish_data_fwd = NULL, model = "IPM", pool_pops = TRUE, 
                    chains = 3, iter = 1000, warmup = 500, seed = 5432,
                    control = list(adapt_delta = 0.95, stepsize = 0.01, max_treedepth = 13))

print(IPM_pp, pars = c("phi","p_HOS","B_rate_all","q","gamma","p","S_tot","R_tot"), include = FALSE)
launch_shinystan(IPM_pp)

# Model with intervention in 1970
IPM_pp_1970 <- salmonIPM(fish_data = fish_data, env_data = pre_post_1970, model = "IPM", pool_pops = TRUE, 
                    chains = 3, iter = 1000, warmup = 500,
                    control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

print(IPM_pp_1970, pars = c("phi","p_HOS","B_rate_all","q","gamma","p","S_tot","R_tot"), include = FALSE)
launch_shinystan(IPM_pp_1970)



#===========================================================================
# FIT PVA MODELS
#===========================================================================

#------------------------------------------
# RR, each population separate
#------------------------------------------

PVA_RR_npp <- salmonIPM(fish_data = fish_data_aug, model = "RR", pool_pops = FALSE, 
                        chains = 3, iter = 1000, warmup = 500,
                        control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

print(PVA_RR_npp, pars = c("phi","R_hat","S_sim"), include = FALSE)
launch_shinystan(PVA_RR_npp)


#------------------------------------------
# RR, populations hierarchical
#------------------------------------------

PVA_RR_pp <- salmonIPM(fish_data = fish_data_aug, model = "RR", pool_pops = TRUE, 
                       chains = 3, iter = 1000, warmup = 500,
                       control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

print(PVA_RR_pp, pars = c("phi","R_hat","S_sim"), include = FALSE)
launch_shinystan(PVA_RR_pp)


#------------------------------------------
# IPM, each population separate
#------------------------------------------

PVA_IPM_npp <- salmonIPM(fish_data = fish_data_aug, model = "IPM", pool_pops = FALSE, 
                         chains = 3, iter = 1000, warmup = 500,
                         control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

print(PVA_IPM_npp, pars = c("phi","p_HOS","B_rate_all","q"), include = FALSE)
launch_shinystan(PVA_IPM_npp)


#------------------------------------------
# IPM, populations hierarchical
#------------------------------------------

PVA_IPM_pp <- salmonIPM(fish_data = fish_data_aug, model = "IPM", pool_pops = TRUE, 
                        chains = 3, iter = 1000, warmup = 500,
                        control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

print(PVA_IPM_pp, pars = c("phi","p_HOS","B_rate_all","q","gamma","p","S_tot","R_tot"), include = FALSE)
launch_shinystan(PVA_IPM_pp)


#===========================================================================
# SIMULATE DATA AND FIT
#===========================================================================

# Create output objects
N_sim <- 30  # number of simulated datasets
N_pop <- length(levels(fish_data$pop))
RR_npp_fit_sim_pars <- list(median = NULL, CI.025 = NULL, CI.975 = NULL)
RR_pp_fit_sim_pars <- list(median = NULL, CI.025 = NULL, CI.975 = NULL)
IPM_npp_fit_sim_pars <- list(median = NULL, CI.025 = NULL, CI.975 = NULL)
IPM_pp_fit_sim_pars <- list(median = NULL, CI.025 = NULL, CI.975 = NULL)
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
  
  IPM_npp_fit_sim_pars[[i]] <- as.data.frame(cbind(matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("a[", 1:N_pop, "]"))),
                                                   matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("Rmax[", 1:N_pop, "]"))),
                                                   matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("rho_proc[", 1:N_pop, "]"))),
                                                   matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("sigma_proc[", 1:N_pop, "]"))),
                                                   matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("sigma_obs[", 1:N_pop, "]")))))
  IPM_pp_fit_sim_pars[[i]] <- as.data.frame(cbind(mu_log_a = rep(NA,N_sim), sigma_log_a = NA,
                                                  matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("a[", 1:N_pop, "]"))),
                                                  mu_log_Rmax = NA, sigma_log_Rmax = NA, 
                                                  matrix(NA, N_sim, N_pop, dimnames = list(NULL, paste0("Rmax[", 1:N_pop, "]"))),
                                                  rho_log_aRmax = NA, rho_log_phi = NA, sigma_log_phi = NA, sigma_proc = NA, sigma_obs = NA,
                                                  matrix(NA, N_sim, 3, dimnames = list(NULL,c("mu_p[1]","mu_p[2]","mu_p[3]"))),
                                                  matrix(NA, N_sim, 2, dimnames = list(NULL,c("sigma_gamma[1]","sigma_gamma[2]"))),
                                                  matrix(NA, N_sim, 2, dimnames = list(NULL,c("sigma_alr_p[1]","sigma_alr_p[2]")))))
}

# Loop over N_sim
for(i in 1:N_sim)
{
  cat("simulated dataset", i, "of", N_sim, "\n")
  
  # Simulate data
  dat <- stan_data(fish_data, model = "IPM")
  sim_fish_data <- IPM_adult_sim(pars = list(mu_log_a = 2, sigma_log_a = 0.5,
                                             mu_log_Rmax = 1, sigma_log_Rmax = 0.8,
                                             rho_log_aRmax = 0.7,
                                             beta_log_phi = 0, rho_log_phi = 0.8, sigma_log_phi = 0.5,
                                             sigma_proc = 0.5, mu_p = c(0.1,0.5,0.4),
                                             sigma_gamma = c(0.1,0.3), L_gamma = diag(2),
                                             sigma_alr_p = c(0.2,0.2), L_alr_p = diag(2),
                                             p_HOS = pmin(fish_data$n_H_obs/(fish_data$n_H_obs + fish_data$n_W_obs), 0.9),
                                             sigma_obs = 0.5),
                                 pop = dat$pop,
                                 year = dat$year,
                                 N_age = 3,
                                 max_age = 5,
                                 A = dat$A,
                                 n_age_tot_obs = rowSums(dat$n_age_obs),
                                 n_HW_tot_obs = fish_data$n_H_obs + fish_data$n_W_obs,
                                 F_rate = dat$F_rate,
                                 B_rate = fish_data$B_take_obs/(fish_data$B_take_obs + dat$S_tot_obs))
  
  sim_fish_data$sim_dat$S_tot_obs[is.na(fish_data$S_tot_obs)] <- NA
  
  # Fit single-pop spawner-recruit models and store estimates
  RR_npp_fit_sim <- salmonIPM(fish_data = sim_fish_data$sim_dat, model = "RR", pool_pops = FALSE,
                             chains = 3, iter = 1000, warmup = 500,
                             control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))
  
  RR_npp_fit_sim <- as.matrix(RR_npp_fit_sim, names(RR_npp_fit_sim_pars$median))
  RR_npp_fit_sim_pars$median[i,] <- apply(RR_npp_fit_sim, 2, median)
  RR_npp_fit_sim_pars$CI.025[i,] <- apply(RR_npp_fit_sim, 2, quantile, 0.025)
  RR_npp_fit_sim_pars$CI.975[i,] <- apply(RR_npp_fit_sim, 2, quantile, 0.975)
  
  # Fit hierarchical spawner-recruit model and store estimates
  RR_pp_fit_sim <- salmonIPM(fish_data = sim_fish_data$sim_dat, model = "RR", 
                          chains = 3, iter = 1000, warmup = 500,
                          control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))
  
  RR_pp_fit_sim <- as.matrix(RR_pp_fit_sim, names(RR_pp_fit_sim_pars$median))
  RR_pp_fit_sim_pars$median[i,] <- apply(RR_pp_fit_sim, 2, median)
  RR_pp_fit_sim_pars$CI.025[i,] <- apply(RR_pp_fit_sim, 2, quantile, 0.025)
  RR_pp_fit_sim_pars$CI.975[i,] <- apply(RR_pp_fit_sim, 2, quantile, 0.975)
  
  # Fit single-pop IPMs and store estimates
  IPM_npp_fit_sim <- salmonIPM(fish_data = sim_fish_data$sim_dat, model = "IPM", pool_pops = FALSE,
                              chains = 3, iter = 2000, warmup = 1000, 
                              control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))
  
  IPM_npp_fit_sim <- as.matrix(IPM_npp_fit_sim, names(IPM_npp_fit_sim_pars$median))
  IPM_npp_fit_sim_pars$median[i,] <- apply(IPM_npp_fit_sim, 2, median)
  IPM_npp_fit_sim_pars$CI.025[i,] <- apply(IPM_npp_fit_sim, 2, quantile, 0.025)
  IPM_npp_fit_sim_pars$CI.975[i,] <- apply(IPM_npp_fit_sim, 2, quantile, 0.975)
  
  # Fit hierarchical IPM and store estimates
  IPM_pp_fit_sim <- salmonIPM(fish_data = sim_fish_data$sim_dat, model = "IPM", 
                           chains = 3, iter = 2000, warmup = 1000, 
                           control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))
  
  IPM_pp_fit_sim <- as.matrix(IPM_pp_fit_sim, names(IPM_pp_fit_sim_pars$median))
  IPM_pp_fit_sim_pars$median[i,] <- apply(IPM_pp_fit_sim, 2, median)
  IPM_pp_fit_sim_pars$CI.025[i,] <- apply(IPM_pp_fit_sim, 2, quantile, 0.025)
  IPM_pp_fit_sim_pars$CI.975[i,] <- apply(IPM_pp_fit_sim, 2, quantile, 0.975)
}


# Plot estimated and true values
dev.new(width = 10, height = 10)
par(mfrow = c(4,4), mar = c(2,3,2,1))
c1 <- "orangered3"
c1 <- col2rgb(c1)
c1 <- rgb(c1[1], c1[2], c1[3], maxColorValue = 255, alpha = 255*0.5)
c2 <- "blue4"
c2 <- col2rgb(c2)
c2 <- rgb(c2[1], c2[2], c2[3], maxColorValue = 255, alpha = 255*0.5)
for(i in names(IPM_pp_fit_sim_pars$median)[1:16])
{
  theta <- unlist(strsplit(strsplit(i, "[", fixed = T)[[1]], "]", fixed = T))
  theta <- sim_fish_data$pars_out[[theta[1]]][ifelse(length(theta)==1, 1, as.numeric(theta[2]))]
  RR_pp_median <- RR_pp_fit_sim_pars$median[,i] - theta
  RR_pp_CI.025 <- RR_pp_fit_sim_pars$CI.025[,i] - theta
  RR_pp_CI.975 <- RR_pp_fit_sim_pars$CI.975[,i] - theta
  IPM_pp_median <- IPM_pp_fit_sim_pars$median[,i] - theta
  IPM_pp_CI.025 <- IPM_pp_fit_sim_pars$CI.025[,i] - theta
  IPM_pp_CI.975 <- IPM_pp_fit_sim_pars$CI.975[,i] - theta
  
  plot(RR_pp_median, IPM_pp_median, pch = "", las = 1, cex.axis = 1.2, cex.lab = 1.5,
       xlim = range(RR_pp_CI.025, RR_pp_CI.975, IPM_pp_CI.025, IPM_pp_CI.975),
       ylim = range(RR_pp_CI.025, RR_pp_CI.975, IPM_pp_CI.025, IPM_pp_CI.975),
       xlab = "", ylab = "", main = i)
  abline(0, 1, col = "darkgray", lwd = 2)
  abline(v = 0, h = 0, col = "darkgray", lwd = 2)
  segments(RR_pp_CI.025, IPM_pp_median, x1 = RR_pp_CI.975, col = c1)
  segments(RR_pp_median, IPM_pp_CI.025, y1 = IPM_pp_CI.975, col = c2)
  points(RR_pp_median, IPM_pp_median, pch = 16, cex = 1.2)
}

rm(list = c("c1","c2","RR_pp_median","RR_pp_CI.025","RR_pp_CI.975",
            "IPM_pp_median","RR_pp_CI.025","RR_pp_CI.975"))


#===========================================================================
# FIGURES
#===========================================================================

#-----------------------------------------------------------------------------------------
# Time series of observed and fitted/predicted total spawners and R/S for 3 populations
#-----------------------------------------------------------------------------------------

dev.new(width=11,height=7)
# png(filename="Fig_1.png", width=11, height=7, units="in", res=200, type="cairo-png")
par(mfcol=c(2,3), mar=c(1,3,4.1,1), oma=c(5.1,2.1,0,0))
pops <- c("Marsh","Catherine","Yankee")
S_tot_IPM <- extract1(PVA_IPM_pp,"S_tot")
S_tot_obs_IPM <- S_tot_IPM * rlnorm(length(S_tot_IPM), 0, extract1(PVA_IPM_pp,"sigma_obs"))
R_tot_IPM <- extract1(PVA_IPM_pp,"R_tot")
RS_IPM <- R_tot_IPM/S_tot_IPM
sd <- stan_data(fish_data_aug, model = "RR")
S_tot_RR <- extract1(PVA_RR_pp,"S_sim")
S_tot_RR[,fish_data_aug$type=="past"] <- NA
RS_RR <- extract1(PVA_RR_pp,"R_sim")
RS_RR[,sd$which_fit] <- extract1(PVA_RR_pp,"R_hat")[,sd$which_fit]
RS_RR <- RS_RR/extract1(PVA_RR_pp,"S_sim")
RS_RR[,fish_data_aug$type=="past" & sd$S_NA] <- NA
Y <- 10
c1 <- "blue4"
c1t <- col2rgb(c1)
c1tt <- c1t
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.3)
c1tt <- rgb(c1tt[1], c1tt[2], c1tt[3], maxColorValue = 255, alpha = 255*0.2)
c2 <- "orangered3"
c2 <- col2rgb(c2)
c2 <- rgb(c2[1], c2[2], c2[3], maxColorValue = 255, alpha = 255*0.6)
for(i in pops)
{
  y1 <- fish_data$year[fish_data$pop==i]
  y2 <- c(y1, max(y1) + 1:Y) 
  plot(y1, fish_data$S_tot_obs[fish_data$pop==i], pch = "",
       xlim = range(fish_data$year[fish_data$pop %in% pops]) + c(0,Y),
       ylim = range(pmax(fish_data$S_tot_obs[fish_data$pop==i], 1),
                    apply(S_tot_obs_IPM[,fish_data_aug$pop %in% pops & fish_data_aug$year <= max(y2)], 2, quantile, c(0.025,0.975)), 
                    na.rm = T), 
       cex.axis = 1.5, cex.main = 2, las = 1, yaxt = "n",
       xlab = "", ylab = "", main = i, log = "y")
  at <- maglab(10^par("usr")[3:4], log = T)
  axis(2, at$labat, cex.axis=1.5, las=1,
       labels = sapply(log10(at$labat), function(i) as.expression(bquote(10^ .(i)))))
  if(i==pops[1]) mtext("Spawners", side = 2, line = 3.5, cex = par("cex")*2)
  lines(y2, apply(S_tot_IPM[,fish_data_aug$pop==i & fish_data_aug$year %in% y2], 2, median), col = c1, lwd = 2)
  polygon(c(y2, rev(y2)), 
          c(apply(S_tot_IPM[,fish_data_aug$pop==i & fish_data_aug$year %in% y2], 2, quantile, 0.025), 
            rev(apply(S_tot_IPM[,fish_data_aug$pop==i & fish_data_aug$year %in% y2], 2, quantile, 0.975))),
          col = c1t, border = NA)
  polygon(c(y2, rev(y2)), 
          c(apply(S_tot_obs_IPM[,fish_data_aug$pop==i & fish_data_aug$year %in% y2], 2, quantile, 0.025), 
            rev(apply(S_tot_obs_IPM[,fish_data_aug$pop==i & fish_data_aug$year %in% y2], 2, quantile, 0.975))),
          col = c1tt, border = NA)
  points(y2, apply(S_tot_RR[,fish_data_aug$pop==i & fish_data_aug$year %in% y2], 2, median, na.rm = T), pch = 16, cex = 1.5, col = c2)
  segments(x0 = y2,
           y0 = apply(S_tot_RR[,fish_data_aug$pop==i & fish_data_aug$year %in% y2], 2, quantile, 0.025, na.rm = T),
           y1 = apply(S_tot_RR[,fish_data_aug$pop==i & fish_data_aug$year %in% y2], 2, quantile, 0.975, na.rm = T),
           col = c2)
  points(y1, fish_data$S_tot_obs[fish_data$pop==i], pch=16, cex = 1.5)
  
  plot(y2, apply(RS_IPM[,fish_data_aug$pop==i & fish_data_aug$year %in% y2], 2, median), pch = "",
       xlim = range(fish_data$year[fish_data$pop %in% pops]) + c(0,Y),
       ylim = range(apply(RS_IPM[,fish_data_aug$pop %in% pops & fish_data_aug$year <= max(y2)], 2, quantile, c(0.025,0.975)),
                    apply(RS_RR[,fish_data_aug$pop %in% pops & fish_data_aug$year <= max(y2)], 2, quantile, c(0.025,0.975), na.rm = T), 
                    na.rm = T), 
       cex.axis = 1.5, las = 1, yaxt = "n",
       xlab = "", ylab = "", log = "y")
  at <- maglab(10^par("usr")[3:4], log = T)
  axis(2, at$labat, cex.axis=1.5, las=1,
       labels = sapply(log10(at$labat), function(i) as.expression(bquote(10^ .(i)))))
  mtext("Year", side = 1, line = 3, cex = par("cex")*2)
  if(i==pops[1]) mtext("Recruits per spawner", side = 2, line = 3.5, cex = par("cex")*2)
  abline(h = 1, lty = 2, lwd = 2, col = "gray")
  lines(y2, apply(RS_IPM[,fish_data_aug$pop==i & fish_data_aug$year %in% y2], 2, median), lwd = 2, col = c1)
  polygon(c(y2, rev(y2)), 
          c(apply(RS_IPM[,fish_data_aug$pop==i & fish_data_aug$year %in% y2], 2, quantile, 0.025), 
            rev(apply(RS_IPM[,fish_data_aug$pop==i & fish_data_aug$year %in% y2], 2, quantile, 0.975))),
          col = c1t, border = NA)
  points(y2, apply(RS_RR[,fish_data_aug$pop==i & fish_data_aug$year %in% y2], 2, median, na.rm = T), pch = 16, cex = 1.5, col = c2)
  segments(x0 = y2, 
           y0 = apply(RS_RR[,fish_data_aug$pop==i & fish_data_aug$year %in% y2], 2, quantile, 0.025, na.rm = T), 
           y1 = apply(RS_RR[,fish_data_aug$pop==i & fish_data_aug$year %in% y2], 2, quantile, 0.975, na.rm = T), 
           col = c2)
}

rm(list = c("pops","S_tot_IPM","S_tot_obs_IPM","R_tot_IPM","RS_IPM","RS_RR","sd","at",
            "c1","c1t","c1tt","c2","Y","y1","y2"))
# dev.off()


#--------------------------------------------------------------------
# Comparison of S-R curves and parameters under RR and IPM models
#--------------------------------------------------------------------

dev.new(width = 11, height = 3.5)
# png(filename="Fig_2.png", width=11, height=3.5, units="in", res=200, type="cairo-png")
par(mfrow = c(1,3), mar = c(5.1,5.1,1,1))
BH <- function(a, Rmax, S) 
{
  a*S/(1 + a*S/Rmax)
}

S <- matrix(seq(0, quantile(fish_data$S_tot_obs/fish_data$A, 0.9, na.rm = T), length = 100),
            nrow = sum(PVA_RR_pp@sim$n_save - PVA_RR_pp@sim$warmup2), ncol = 100, byrow = T)

# S-R curves
mu_log_a <- as.vector(extract1(PVA_RR_pp,"mu_log_a"))
mu_log_Rmax <- as.vector(extract1(PVA_RR_pp,"mu_log_Rmax"))
R_ESU_RR <- BH(a = exp(mu_log_a), Rmax = exp(mu_log_Rmax), S = S)
mu_log_a <- as.vector(extract1(PVA_IPM_pp,"mu_log_a"))
mu_log_Rmax <- as.vector(extract1(PVA_IPM_pp,"mu_log_Rmax"))
R_ESU_IPM <- BH(a = exp(mu_log_a), Rmax = exp(mu_log_Rmax), S = S)

bb <- "orangered3"
plot(S[1,], apply(R_ESU_RR, 2, median), type = "l", lwd=3, col = bb, las = 1,
     cex.lab = 2, cex.axis = 1.5, xaxs = "i", yaxs = "i",
     ylim = range(0, apply(R_ESU_RR, 2, quantile, 0.975), apply(R_ESU_IPM, 2, quantile, 0.975)),
     xlab=bquote("Spawner density (ha"^-1*")"), ylab=bquote("Recruit density (ha"^-1*")"))
bb <- col2rgb(bb)
bb <- rgb(bb[1], bb[2], bb[3], maxColorValue = 255, alpha = 255*0.15)
polygon(c(S[1,], rev(S[1,])), c(apply(R_ESU_RR, 2, quantile, 0.025), rev(apply(R_ESU_RR, 2, quantile, 0.975))), col = bb, border = NA)

bb <- "blue4"
lines(S[1,], apply(R_ESU_IPM, 2, median), type = "l", lwd=3, col = bb)
bb <- col2rgb(bb)
bb <- rgb(bb[1], bb[2], bb[3], maxColorValue = 255, alpha = 255*0.3)
polygon(c(S[1,], rev(S[1,])), c(apply(R_ESU_IPM, 2, quantile, 0.025), rev(apply(R_ESU_IPM, 2, quantile, 0.975))), col = bb, border = NA)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "A", cex = 2)

# Posterior densities of log(a)
dd_IPM_ESU <- density(extract1(PVA_IPM_pp,"mu_log_a"))
dd_IPM_pop <- vector("list", length(levels(fish_data$pop)))
for(i in 1:length(dd_IPM_pop))
  dd_IPM_pop[[i]] <- density(log(extract1(PVA_IPM_pp,"a")[,i]))
dd_RR_ESU <- density(extract1(PVA_RR_pp,"mu_log_a"))
dd_RR_pop <- vector("list", length(levels(fish_data$pop)))
for(i in 1:length(dd_RR_pop))
  dd_RR_pop[[i]] <- density(log(extract1(PVA_RR_pp,"a")[,i]))

bb <- "blue4"
plot(dd_IPM_ESU$x, dd_IPM_ESU$y, type = "l", lwd = 3, col = bb, las = 1, cex.lab = 2, cex.axis = 1.5,
     xlab = bquote(log(alpha)), ylab = "Probability density", xaxs = "i",
     xlim = range(c(dd_IPM_ESU$x, dd_RR_ESU$x, sapply(c(dd_RR_pop, dd_IPM_pop), function(m) m$x))),
     ylim = range(c(dd_IPM_ESU$y, dd_RR_ESU$y, sapply(c(dd_RR_pop, dd_IPM_pop), function(m) m$y))))
bb <- col2rgb(bb)
bb <- rgb(bb[1], bb[2], bb[3], maxColorValue = 255, alpha = 255*0.3)
for(i in 1:length(dd_IPM_pop))
  lines(dd_IPM_pop[[i]]$x, dd_IPM_pop[[i]]$y, col = bb)

bb <- "orangered3"
lines(dd_RR_ESU$x, dd_RR_ESU$y, lwd = 3, col = bb)
bb <- col2rgb(bb)
bb <- rgb(bb[1], bb[2], bb[3], maxColorValue = 255, alpha = 255*0.3)
for(i in 1:length(dd_RR_pop))
  lines(dd_RR_pop[[i]]$x, dd_RR_pop[[i]]$y, col = bb)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "B", cex = 2)

# Posterior densities of log(Rmax)
dd_IPM_ESU <- density(extract1(PVA_IPM_pp,"mu_log_Rmax"))
dd_IPM_pop <- vector("list", length(levels(fish_data$pop)))
for(i in 1:length(dd_IPM_pop))
  dd_IPM_pop[[i]] <- density(log(extract1(PVA_IPM_pp,"Rmax")[,i]))
dd_RR_ESU <- density(extract1(PVA_RR_pp,"mu_log_Rmax"))
dd_RR_pop <- vector("list", length(levels(fish_data$pop)))
for(i in 1:length(dd_RR_pop))
  dd_RR_pop[[i]] <- density(log(extract1(PVA_RR_pp,"Rmax")[,i]))

bb <- "blue4"
plot(dd_IPM_ESU$x, dd_IPM_ESU$y, type = "l", lwd = 3, col = bb, las = 1, cex.lab = 2, cex.axis = 1.5,
     xlab = bquote(log(italic(R)[max]*" [ha"^-1*"]")), ylab = "Probability density", xaxs = "i",
     xlim = range(c(dd_IPM_ESU$x, dd_RR_ESU$x, sapply(c(dd_RR_pop, dd_IPM_pop), function(m) m$x))),
     ylim = range(c(dd_IPM_ESU$y, dd_RR_ESU$y, sapply(c(dd_RR_pop, dd_IPM_pop), function(m) m$y))))
bb <- col2rgb(bb)
bb <- rgb(bb[1], bb[2], bb[3], maxColorValue = 255, alpha = 255*0.3)
for(i in 1:length(dd_IPM_pop))
  lines(dd_IPM_pop[[i]]$x, dd_IPM_pop[[i]]$y, col = bb)

bb <- "orangered3"
lines(dd_RR_ESU$x, dd_RR_ESU$y, lwd = 3, col = bb)
bb <- col2rgb(bb)
bb <- rgb(bb[1], bb[2], bb[3], maxColorValue = 255, alpha = 255*0.3)
for(i in 1:length(dd_RR_pop))
  lines(dd_RR_pop[[i]]$x, dd_RR_pop[[i]]$y, col = bb)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "C", cex = 2)
legend("topright", c("IPM","RR"), col = c("blue4","orangered3"), lwd = 3, cex = 1.5)

rm(list=c("mu_log_a","mu_log_Rmax","S","R_ESU_RR","R_ESU_IPM","BH",
          "bb","dd_IPM_ESU","dd_RR_ESU","dd_IPM_pop","dd_RR_pop"))
# dev.off()


#-------------------------------------------------------------------------
# Pop-specific S-R curves with and without obs error
#-------------------------------------------------------------------------

dev.new(width = 7, height = 7)
# png(filename="Fig_2.5(3).png", width=7, height=7, units="in", res=200, type="cairo-png")
BH <- function(a, Rmax, S, A) 
{
  a*S/(A + a*S/Rmax)
}

pop <- "Sulphur"

yy <- stan_data(fish_data_aug, model = "RR")$year
AA <- unique(fish_data$A[fish_data$pop==pop])
S_tot_RR <- fish_data$S_tot_obs[fish_data$pop==pop]
R_adj_RR <- colMeans(sweep(1/extract1(PVA_RR_pp,"phi")[,yy], 2, run_recon(fish_data_aug)$R, "*"))[fish_data_aug$pop==pop & fish_data_aug$type=="past"]
S_tot_IPM <- extract1(PVA_IPM_pp,"S_tot")[,fish_data_aug$pop==pop & fish_data_aug$type=="past"]
R_adj_IPM <- (extract1(PVA_IPM_pp,"R_tot")/extract1(PVA_IPM_pp,"phi")[,yy])[,fish_data_aug$pop==pop & fish_data_aug$type=="past"]

S <- matrix(seq(0, max(S_tot_RR, apply(S_tot_IPM, 2, median), na.rm = T)*1.02, length = 500),
            nrow = sum(PVA_RR_pp@sim$n_save - PVA_RR_pp@sim$warmup2), ncol = 500, byrow = T)
a <- as.vector(extract1(PVA_RR_pp,"a")[,which(levels(fish_data$pop)==pop)])
Rmax <- as.vector(extract1(PVA_RR_pp,"Rmax")[,which(levels(fish_data$pop)==pop)])
R_RR <- BH(a = a, Rmax = Rmax, S = S, AA) * AA
a <- as.vector(extract1(PVA_IPM_pp,"a")[,which(levels(fish_data$pop)==pop)])
Rmax <- as.vector(extract1(PVA_IPM_pp,"Rmax")[,which(levels(fish_data$pop)==pop)])
R_IPM <- BH(a = a, Rmax = Rmax, S = S, AA) * AA

# layer 1
c1 <- "orangered3"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.15)
plot(S[1,], apply(R_RR, 2, median), type = "l", lwd=3, col = c1, las = 1,
     cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5, xaxs = "i", yaxs = "i",
     xlim = c(-3, max(S)),
     ylim = range(0, R_adj_RR, apply(R_adj_IPM, 2, quantile, 0.975), apply(R_RR, 2, quantile, 0.975), apply(R_IPM, 2, quantile, 0.975), na.rm = T)*1.05,
     xlab = "Spawners", ylab = "Adjusted recruits", main = pop)
polygon(c(S[1,], rev(S[1,])), c(apply(R_RR, 2, quantile, 0.025), rev(apply(R_RR, 2, quantile, 0.975))), col = c1t, border = NA)
points(S_tot_RR, R_adj_RR, pch = 1)

# layer 2
c2 <- "blue4"
c2t <- col2rgb(c2)
c2tt <- c2t
c2t <- rgb(c2t[1], c2t[2], c2t[3], maxColorValue = 255, alpha = 255*0.4)
c2tt <- rgb(c2tt[1], c2tt[2], c2tt[3], maxColorValue = 255, alpha = 255*0.3)
points(apply(S_tot_IPM, 2, median), apply(R_adj_IPM, 2, median), pch = 16, col = c2)
arrows(S_tot_RR, R_adj_RR, apply(S_tot_IPM, 2, median), apply(R_adj_IPM, 2, median), col = c2t, length = 0.1)
segments(x0 = apply(S_tot_IPM, 2, quantile, 0.025), y0 = apply(R_adj_IPM, 2, median), x1 = apply(S_tot_IPM, 2, quantile, 0.975), col = c2t)
segments(x0 = apply(S_tot_IPM, 2, median), y0 = apply(R_adj_IPM, 2, quantile, 0.025), y1 = apply(R_adj_IPM, 2, quantile, 0.975), col = c2t)

# layer 3
lines(S[1,], apply(R_IPM, 2, median), lwd = 3, col = c2)
polygon(c(S[1,], rev(S[1,])), c(apply(R_IPM, 2, quantile, 0.025), rev(apply(R_IPM, 2, quantile, 0.975))), col = c2tt, border = NA)

rm(list = c("yy","AA","BH","pop","S","a","Rmax","R_RR","R_IPM","S_tot_RR","S_tot_IPM",
            "R_adj_RR","R_adj_IPM","c1","c1t","c2","c2t","c2tt"))

# dev.off()


#------------------------------------------------------------------
# Probability of quasi-extinction by population under RR and IPM
#------------------------------------------------------------------

dev.new(height = 7, width = 7)
# png(filename="Fig_3.png", width=7, height=7, units="in", res=200, type="cairo-png")
par(oma = c(0,7,0,0))
qet <- 50     # set quasi-extinction threshold (4-yr moving average)
pop <- fish_data_aug$pop[fish_data_aug$type=="future"]
S_tot_RR <- t(extract1(PVA_RR_pp, "S_sim")[,fish_data_aug$type=="future"])
pqe_RR <- aggregate(S_tot_RR, list(pop = pop), function(x) any(rollmean(x, 4) < qet))
S_tot_IPM <- t(extract1(PVA_IPM_pp, "S_tot")[,fish_data_aug$type=="future"])
pqe_IPM <- aggregate(S_tot_IPM, list(pop = pop), function(x) any(rollmean(x, 4) < qet))
pqe <- data.frame(pop = pqe_IPM[,1], pqe_RR = rowMeans(pqe_RR[,-1]), pqe_IPM = rowMeans(pqe_IPM[,-1]))
pqe <- pqe[order(pqe$pqe_IPM),]

barplot(t(pqe[,-1]), names.arg = pqe$pop, horiz = TRUE, beside = TRUE,
        col = c("orangered3","blue4"), las = 1, cex.axis = 1, cex.lab = 1.2, cex.names = 0.9, cex.main = 1.2,
        xlab = "Probability of quasi-extinction")
# main = paste("50-year quasi-extinction risk \n QET = ", qet, " spawners (4-year moving average)", sep = ""))
legend("bottomright", inset = 0.1, c("IPM","RR"), pch = 15, pt.cex = 2, cex = 1, col = c("blue4","orangered3"))

rm(list=c("qet","pop","S_tot_RR","S_tot_IPM","pqe_RR","pqe_IPM","pqe"))
# dev.off()


#-------------------------------------------------------------------------
# ESU-level probabilities of quasi-extinction as a function of QET
#-------------------------------------------------------------------------

dev.new(height = 7, width = 7)
# png(filename="Fig_4.png", width=7, height=7, units="in", res=200, type="cairo-png")
par(mar = c(5.1,5.1,1,1))
qet <- 0:50 
pop <- fish_data_aug$pop[fish_data_aug$year > max(fish_data$year)]
year <- fish_data_aug$year[fish_data_aug$year > max(fish_data$year)]
S_tot_RR <- extract1(PVA_RR_pp, "S_sim")[,fish_data_aug$year > max(fish_data$year)]
S_tot_RR <- tapply(as.vector(S_tot_RR), 
                   list(iter = rep(1:nrow(S_tot_RR), ncol(S_tot_RR)), 
                        year = rep(year, each = nrow(S_tot_RR)), 
                        pop = rep(pop, each = nrow(S_tot_RR))), 
                   identity)
S_tot_RR <- apply(S_tot_RR, c(1,3), function(x) rollmean(x, 4))
pqe_RR <- sapply(qet, function(qq) mean(apply(apply(S_tot_RR, c(2,3), function(x) any(x < qq)), 1, any)))

S_tot_IPM <- extract1(PVA_IPM_pp, "S_tot")[,fish_data_aug$year > max(fish_data$year)]
S_tot_IPM <- tapply(as.vector(S_tot_IPM), 
                    list(iter = rep(1:nrow(S_tot_IPM), ncol(S_tot_IPM)), 
                         year = rep(year, each = nrow(S_tot_IPM)), 
                         pop = rep(pop, each = nrow(S_tot_IPM))), 
                    identity)
S_tot_IPM <- apply(S_tot_IPM, c(1,3), function(x) rollmean(x, 4))
pqe_IPM <- sapply(qet, function(qq) mean(apply(apply(S_tot_IPM, c(2,3), function(x) any(x < qq)), 1, any)))

plot(qet, pqe_RR, type = "l", lwd = 3, col = "orangered3",
     ylim = range(pqe_RR, pqe_IPM), xaxs = "i", yaxs = "i", cex.axis = 1.2, cex.lab = 1.5, las = 1,
     xlab = "Quasi-extinction threshold", ylab = bquote("Probability of ">="1 quasi-extinction"))
lines(qet, pqe_IPM, lwd = 3, col = "blue4")
legend("topleft", c("IPM","RR"), col = c("blue4","orangered3"), lwd = 3, cex = 1.5)

rm(list = c("pop","S_tot_RR","S_tot_IPM","qet","year","pqe_RR","pqe_IPM"))
# dev.off()


#------------------------------------------------------------------------------
# CDF of max sustainable harvest rate, at ESU and pop levels, under RR and IPM
# Umax = 1 - 1/a
#------------------------------------------------------------------------------

dev.new(width = 7, height = 8)
# png(filename="Fig_5.png", width=7, height=8, units="in", res=200, type="cairo-png")
layout(matrix(c(2,1)), heights = c(1.5,7))

mu_log_a_RR <- extract1(PVA_RR_pp,"mu_log_a")
Umax_ESU_RR <- 1 - exp(-mu_log_a_RR)
a_RR <- extract1(PVA_RR_pp,"a")
Umax_pop_RR <- 1 - 1/a_RR

mu_log_a_IPM <- extract1(PVA_IPM_pp,"mu_log_a")
Umax_ESU_IPM <- 1 - exp(-mu_log_a_IPM)
a_IPM <- extract1(PVA_IPM_pp,"a")
Umax_pop_IPM <- 1 - 1/a_IPM

M <- length(mu_log_a_RR)

c1 <- "orangered3"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.6)
c2 <- "blue4"
c2t <- col2rgb(c2)
c2t <- rgb(c2t[1], c2t[2], c2t[3], maxColorValue = 255, alpha = 255*0.7)

par(mar = c(5.1,4.1,0,2.1))

plot(sort(Umax_ESU_RR), (1:M)/M, type = "l", lwd = 4, col = c1,
     xlim = c(0, 1), ylim = c(0,1),
     xaxs = "i", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     xlab = "Mortality rate", ylab = "Probability of decline")
lines(sort(Umax_ESU_IPM), (1:M)/M, lwd = 4, col = c2)
for(i in 1:ncol(Umax_pop_IPM))
{
  lines(sort(Umax_pop_RR[,i]), (1:M)/M, lwd = 1, col = c1t)
  lines(sort(Umax_pop_IPM[,i]), (1:M)/M, lwd = 1, col = c2t)
}
legend("topleft", c("IPM (populations)","IPM (metapopulation)","RR (populations)","RR (metapopulation)"), 
       col = c("blue4","blue4","orangered3","orangered3"), lwd = c(1,4,1,4), cex = 1.2)

par(mar = c(0.1,3.1,0.5,1.5))
bins <- seq(0, 1, 0.05)
F1 <- hist(fish_data$F_rate[fish_data$year < 1980], breaks = bins, plot = F)$counts
F2 <- hist(fish_data$F_rate[fish_data$year >= 1980], breaks = bins, plot = F)$counts
barplot(rbind(F1,F2), names = bins[-1], space = 0, yaxs = "i", xaxt = "n", yaxt = "n", 
        xlab = "", ylab = "", main = "", col = c("black","gray"), border = "white")
text(2, max(F2)*0.8, "post-1980 harvest", pos = 4, cex = 1.2, col = "darkgray")
text(8, max(F1)*1.5, "pre-1980 harvest", pos = 4, cex = 1.2)

rm(list = c("mu_log_a_RR","mu_log_a_IPM","Umax_ESU_RR","Umax_ESU_IPM",
            "a_RR","Umax_pop_RR","a_IPM","Umax_pop_IPM","c1","c1t","c2","c2t","bins","F1","F2"))

# dev.off()



#===========================================================================
# MORE FIGURES: SPARE PARTS
#===========================================================================

#--------------------------------------------------------------------------
# Generic example of spawner-recruit regression
#--------------------------------------------------------------------------

dev.new(width = 7, height = 7)
par(mar = c(3,3,1,1))
png(filename="BH_fig.png", width=7, height=7, units="in", res=200, type="cairo-png")
BH <- function(a, b, S)
{
  a*S/(1 + b*S)
}
a <- 1
b <- 1
sigma <- 0.3
S <- seq(0,2,length = 20)
R <- BH(a, b, S)*rlnorm(length(S), -sigma^2/2, sigma)

plot(S, R, xaxt = "n", yaxt = "n", xlab = "", ylab = "", cex.lab = 2, pch = 16, cex = 2)
mtext("Spawners", 1, line = 1, cex = 2)
mtext("Recruits", 2, line = 1, cex = 2)
curve(BH(a, b, x), add = T, col = "blue", lwd = 3)
box()
rm(list = c("BH","a","b","sigma","S","R"))
dev.off()


#--------------------------------------------------------------------------
# Hierarchical S-R curve illustration
#--------------------------------------------------------------------------

png(filename="BH_hierarchical_fig1.png", width=14, height=7, units="in", res=200, type="cairo-png")
par(mfrow = c(1,2))
cols = c("blue", "green", "orange", "red")
curve(dlnorm(x, 0, 0.5), from = 0, to = exp(1), lwd = 5, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
a <- rlnorm(4, 0, 0.5)
segments(x0 = a, y0 = rep(0,4), y1 = dlnorm(a, 0, 0.5), lwd = 5, col = cols) 
mtext(bquote(alpha), 1, line = 3, cex = 4)
mtext("Probability", 2, line = 1, cex = 4)
curve(dlnorm(x, 0, 1), from = 0, to = exp(2), lwd = 5, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
Rmax <- rlnorm(4, 0, 1)
segments(x0 = Rmax, y0 = rep(0,4), y1 = dlnorm(Rmax, 0, 1), lwd = 5, col = cols) 
mtext(bquote(italic(R)[max]), 1, line = 3, cex = 4)
rm(list = c("cols","a","Rmax"))
dev.off()

dev.new(width = 16, height = 4)
png(filename="BH_hierarchical_fig2.png", width=16, height=4, units="in", res=200, type="cairo-png")
par(mfrow=c(1,4))
BH <- function(a, b, S)
{
  a*S/(1 + b*S)
}
a <- rlnorm(4, 0, 0.5)
b <- rlnorm(4, 0, 0.5)
sigma <- 0.3
S <- seq(0,2,length = 20)
R <- matrix(NA, length(S), 4)
for(i in 1:4)
  R[,i] <- BH(a[i], b[i], S)*rlnorm(length(S), -sigma^2/2, sigma)
cols = c("blue", "green", "orange", "red")

for(i in 1:4)
{
  curve(BH(a[i], b[i], x), from = 0, to = max(S), lwd = 3, col = cols[i], ylim = c(0,max(R[,i])),
        xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = paste("Population", i), cex.main = 2.5)
  points(S, R[,i], pch = 16, cex = 2, col = cols[i])
  mtext("Spawners", 1, line = 1.5, cex = 2)
  mtext("Recruits", 2, line = 1, cex = 2)
}
dev.off()

dev.new(width = 7, height = 4)
png(filename="BH_hierarchical_fig3.png", width=7, height=4, units="in", res=200, type="cairo-png")

rho <- 0.7
sigma_log_phi <- 0.3
tt <- 1:30
log_phi <- rep(NA, length(tt))
log_phi[1] <- rnorm(1,0,sigma_log_phi)
for(i in 2:length(tt))
  log_phi[i] <- rnorm(1, log_phi[i-1]*rho, sigma_log_phi)
proc <- matrix(NA, length(tt), 4)
for(i in 1:4)
  proc[,i] <- rnorm(length(tt), log_phi, sigma)

plot(tt, log_phi, type = "l", lwd = 5, xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", ylim = range(proc))
for(i in 1:4)
  lines(tt, proc[,i], lwd = 2, col = cols[i])
lines(tt, log_phi, lwd = 5)
mtext("Brood year", 1, line = 1.5, cex = 2)
mtext("Process error", 2, line = 1, cex = 2)

dev.off()

rm(list = c("BH","a","b","sigma","S","R","rho","sigma_log_phi","tt","log_phi","proc"))


#--------------------------------------------------------------------------
# Spawner-recruit curve variance decomposition under RR and IPM models
#--------------------------------------------------------------------------

dev.new(width=15,height=10)
par(mfcol = c(2,3), mar=c(5.1,5.1,1,1))
BH <- function(a, b, S)
{
  a*S/(1 + b*S)
}

# RR
mu_log_a <- as.vector(extract1(PVA_RR_pp,"mu_log_a"))
sigma_log_a <- as.vector(extract1(PVA_RR_pp,"sigma_log_a"))
mu_log_b <- as.vector(extract1(PVA_RR_pp,"mu_log_b"))
sigma_log_b <- as.vector(extract1(PVA_RR_pp,"sigma_log_b"))
sigma_log_phi <- as.vector(extract1(PVA_RR_pp,"sigma_log_phi"))
sigma <- as.vector(extract1(PVA_RR_pp,"sigma"))
S <- matrix(seq(1, 500, length = 100),
            nrow=sum(PVA_RR_pp@sim$n_save - PVA_RR_pp@sim$warmup2)*50, ncol=100, byrow=T)
R_ESU <- BH(a = exp(mu_log_a), b = exp(mu_log_b), S = S)
R_pop <- BH(a = rlnorm(nrow(S), mu_log_a, sigma_log_a), b = rlnorm(nrow(S), mu_log_b, sigma_log_b), S = S)
R_year <- R_pop*rlnorm(nrow(S), 0, sigma_log_phi)
R_resid <- R_year*rlnorm(nrow(S), 0, sigma)
bb <- "orangered3"
plot(S[1,], apply(R_resid, 2, median), type = "l", lwd=3, col = bb, las = 1,
     cex.lab = 1.5, cex.axis = 1.2, xaxs = "i",
     ylim = c(1, max(apply(R_resid, 2, quantile, 0.975))),
     xlab = "Spawners", ylab = "Recruits", log = "y")
bb <- col2rgb(bb)
bb <- rgb(bb[1], bb[2], bb[3], maxColorValue = 255, alpha = 255*0.15)
polygon(c(S[1,], rev(S[1,])),
        c(apply(R_ESU, 2, quantile, 0.025), rev(apply(R_ESU, 2, quantile, 0.975))), col = bb, border = NA)
polygon(c(S[1,], rev(S[1,])),
        c(apply(R_pop, 2, quantile, 0.025), rev(apply(R_pop, 2, quantile, 0.975))), col = bb, border = NA)
polygon(c(S[1,], rev(S[1,])),
        c(apply(R_year, 2, quantile, 0.025), rev(apply(R_year, 2, quantile, 0.975))), col = bb, border = NA)
polygon(c(S[1,], rev(S[1,])),
        c(apply(R_resid, 2, quantile, 0.025), rev(apply(R_resid, 2, quantile, 0.975))), col = bb, border = NA)

# IPM
mu_log_a <- as.vector(extract1(PVA_IPM_pp,"mu_log_a"))
sigma_log_a <- as.vector(extract1(PVA_IPM_pp,"sigma_log_a"))
mu_log_b <- as.vector(extract1(PVA_IPM_pp,"mu_log_b"))
sigma_log_b <- as.vector(extract1(PVA_IPM_pp,"sigma_log_b"))
sigma_log_phi <- as.vector(extract1(PVA_IPM_pp,"sigma_log_phi"))
mu_sigma_proc <- as.vector(extract1(PVA_IPM_pp,"mu_sigma_proc"))
sigma_obs <- as.vector(extract1(PVA_IPM_pp,"sigma_obs"))
R_ESU <- BH(a = exp(mu_log_a), b = exp(mu_log_b), S = S)
R_pop <- BH(a = rlnorm(nrow(S), mu_log_a, sigma_log_a), b = rlnorm(nrow(S), mu_log_b, sigma_log_b), S = S)
R_year <- R_pop*rlnorm(nrow(S), 0, sigma_log_phi)
R_proc <- R_year*rlnorm(nrow(S), 0, mu_sigma_proc)
R_obs <- R_proc*rlnorm(nrow(S), 0, sigma_obs)
bb <- "blue4"
plot(S[1,], apply(R_obs, 2, median), type = "l", lwd=3, col = bb, las = 1,
     cex.lab = 1.5, cex.axis = 1.2, xaxs = "i",
     ylim = c(1, max(apply(R_obs, 2, quantile, 0.975))),
     xlab = "Spawners", ylab = "Recruits", log = "y")
bb <- col2rgb(bb)
bb <- rgb(bb[1], bb[2], bb[3], maxColorValue = 255, alpha = 255*0.15)
polygon(c(S[1,], rev(S[1,])),
        c(apply(R_ESU, 2, quantile, 0.025), rev(apply(R_ESU, 2, quantile, 0.975))), col = bb, border = NA)
polygon(c(S[1,], rev(S[1,])),
        c(apply(R_pop, 2, quantile, 0.025), rev(apply(R_pop, 2, quantile, 0.975))), col = bb, border = NA)
polygon(c(S[1,], rev(S[1,])),
        c(apply(R_year, 2, quantile, 0.025), rev(apply(R_year, 2, quantile, 0.975))), col = bb, border = NA)
polygon(c(S[1,], rev(S[1,])),
        c(apply(R_proc, 2, quantile, 0.025), rev(apply(R_proc, 2, quantile, 0.975))), col = bb, border = NA)
polygon(c(S[1,], rev(S[1,])),
        c(apply(R_obs, 2, quantile, 0.025), rev(apply(R_obs, 2, quantile, 0.975))), col = bb, border = NA)

rm(list=c("mu_log_a","sigma_log_a","mu_log_b","sigma_log_b","sigma_log_phi","sigma","BH",
          "mu_sigma_proc","sigma_obs","S","R_ESU","R_pop","R_year","R_proc","R_obs","R_resid","bb"))


#--------------------------------------------------------------------------------
# Time series of observed and fitted total spawners for each pop
# under hierarchical IPM
#--------------------------------------------------------------------------------

dev.new(width=16,height=10)
# png(filename="S_tot_fit_IPM.png", width=16*0.9, height=10*0.9, units="in", res=200, type="cairo-png")
par(mfrow=c(5,6), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))

S_tot_IPM <- extract1(IPM_pp,"S_tot")
S_tot_obs_IPM <- S_tot_IPM * rlnorm(length(S_tot_IPM), 0, extract1(IPM_pp,"sigma_obs"))
init_NA <- fish_data$year
for(i in levels(fish_data$code))
  init_NA[fish_data$code==i] <- init_NA[fish_data$code==i] - min(init_NA[fish_data$code==i]) + 1
init_NA <- is.na(fish_data$S_tot_obs) & init_NA < 5

c1 <- "blue4"
c1t <- col2rgb(c1)
c1tt <- c1t
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.3)
c1tt <- rgb(c1tt[1], c1tt[2], c1tt[3], maxColorValue = 255, alpha = 255*0.2)

for(i in levels(fish_data$code))
{
  y1 <- fish_data$year[fish_data$code==i]
  plot(y1, fish_data$S_tot_obs[fish_data$code==i], pch = "",
       xlim = range(fish_data$year),
       ylim = range(pmax(fish_data$S_tot_obs[fish_data$code==i], 1),
                    apply(S_tot_obs_IPM[,fish_data$code==i & !init_NA], 2, quantile, c(0.025,0.975)), na.rm = T), 
       cex.axis = 1.2, las = 1, yaxs = "i", yaxt = "n", xlab = "", ylab = "", log = "y")
  at <- maglab(10^par("usr")[3:4], log = T)
  axis(2, at$labat, cex.axis = 1.2, las = 1,
       labels = sapply(log10(at$labat), function(i) as.expression(bquote(10^ .(i)))))
  mtext(i, side = 3, line = 0.5, cex = par("cex")*1.5)
  if(i %in% levels(fish_data$code)[seq(1,29,6)]) mtext("Spawners", side = 2, line = 3.5, cex = par("cex")*1.5)
  if(i %in% levels(fish_data$code)[25:29]) mtext("Year", side = 1, line = 3, cex = par("cex")*1.5)
  lines(y1, apply(S_tot_IPM[,fish_data$code==i], 2, median), col = c1, lwd = 2)
  polygon(c(y1, rev(y1)), 
          c(apply(S_tot_IPM[,fish_data$code==i], 2, quantile, 0.025), 
            rev(apply(S_tot_IPM[,fish_data$code==i], 2, quantile, 0.975))),
          col = c1t, border = NA)
  polygon(c(y1, rev(y1)), 
          c(apply(S_tot_obs_IPM[,fish_data$code==i], 2, quantile, 0.025), 
            rev(apply(S_tot_obs_IPM[,fish_data$code==i], 2, quantile, 0.975))),
          col = c1tt, border = NA)
  points(y1, fish_data$S_tot_obs[fish_data$code==i], pch=16, cex = 1)
}

rm(list = c("S_tot_IPM","S_tot_obs_IPM","at","c1","c1t","c1tt","y1","init_NA"))
# dev.off()


# #-----------------------------------------------------------------------------------------
# # Time series of observed and fitted or predicted recruits per spawner for each pop
# #-----------------------------------------------------------------------------------------
# 
# dev.new(width=16,height=10)
# par(mfrow=c(4,6), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))
# S_tot <- extract1(PVA_IPM_pp,"S_tot")
# R_tot <- extract1(PVA_IPM_pp,"R_tot")
# RS <- R_tot/S_tot
# rr <- run_recon(fish_data_aug)
# c1 <- "darkgray"
# c1t <- col2rgb(c1)
# c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.6)
# c2 <- "blue"
# c2t <- col2rgb(c2)
# c2t <- rgb(c2t[1], c2t[2], c2t[3], maxColorValue = 255, alpha = 255*0.4)
# iters <- 1:5
# for(i in levels(fish_data_aug$pop))
# {
#   RS_obs_i <- (rr$R/rr$S)[rr$pop==i]
#   RS_i <- RS[,fish_data_aug$pop==i]
#   RS_fit_i <- RS_i[,fish_data_aug$type[fish_data_aug$pop==i]=="past"]
#   RS_fore_i <- cbind(RS_fit_i[,ncol(RS_fit_i)],
#                      RS_i[,fish_data_aug$type[fish_data_aug$pop==i]=="future"])
#   year_fit_i <- fish_data_aug$year[fish_data_aug$pop==i & fish_data_aug$type=="past"]
#   year_fore_i <- c(year_fit_i[length(year_fit_i)], 
#                    fish_data_aug$year[fish_data_aug$pop==i & fish_data_aug$type=="future"])
#   
#   plot(fish_data_aug$year[fish_data_aug$pop==i], RS_obs_i, pch = "", cex=1.2, cex.axis=1.2, las=1,
#        # ylim=c(min(apply(RS_i, 2, quantile, 0.025)),max(apply(RS_i, 2, quantile, 0.975))), 
#        xlim = range(fish_data_aug$year),
#        ylim = range(c(range(RS_i[iters,]),
#                       c(min(apply(RS_fit_i, 2, quantile, 0.025)),max(apply(RS_fit_i, 2, quantile, 0.975))),
#                       range(replace(RS_obs_i, RS_obs_i==0 | is.infinite(RS_obs_i), NA), na.rm = T))),
#        xlab="", ylab="", yaxt = "n", log = "y")
#   at <- maglab(10^par("usr")[3:4], log = T)
#   axis(2, at$labat, cex.axis=1.2, las=1,
#        labels = at$labat)
#   lines(year_fit_i,apply(RS_fit_i, 2, quantile, 0.5), col=c1, lwd=2)
#   polygon(c(year_fit_i, rev(year_fit_i)),
#           c(apply(RS_fit_i, 2, quantile, 0.025), rev(apply(RS_fit_i, 2, quantile, 0.975))),
#           col = c1t, border = NA)
#   # lines(year_fore_i,apply(RS_fore_i, 2, quantile, 0.5), col=c2, lwd=2)
#   # polygon(c(year_fore_i, rev(year_fore_i)),
#   #         c(apply(RS_fore_i, 2, quantile, 0.025), rev(apply(RS_fore_i, 2, quantile, 0.975))),
#   #         col = c2t, border = NA)
#   for(j in iters)
#   {
#     lines(year_fit_i, RS_fit_i[j,], col = c1)
#     lines(year_fore_i, RS_fore_i[j,], col = c2)
#   }
#   points(fish_data_aug$year[fish_data_aug$pop==i], RS_obs_i, type = "b", pch = 16)
#   mtext(i, side=3, line=0.5, cex=1.2)
# }
# mtext("Year", outer = T, side=1, line=2, cex=1.2)
# mtext("Recruits / spawner", outer = T, side=2, line=1.1, cex=1.2)
# rm(list=c("rr","S_tot","RS_obs_i","RS","RS_i","RS_fit_i","RS_fore_i",
#           "year_fit_i","year_fore_i","c1","c1t","c2","c2t","at","iters"))


#-----------------------------------------------------------------------------------------
# Cross-correlation among pops in 
# (1) Recruitment process error residuals log(R_tot) - log(R_tot_hat)
# (2) Spawner observation error residuals log(S_tot_obs) - log(S_tot)
# Correlation matrices are grouped by agglomerative clustering (UPGMA),
# with the largest number of clusters such that each cluster contains >= 3 pops
#-----------------------------------------------------------------------------------------

c1 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
        "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))(200)

# Recruitment residuals
dev.new(width = 10, height = 10, mar = c(1,1,1,1))
# png(filename="log_R_tot_z_corr.png", width=10*0.9, height=10*0.9, units="in", res=200, type="cairo-png")
cor_log_R_tot_z <- cor(tapply(stan_mean(IPM_pp,"log_R_tot_z"), list(fish_data$year, fish_data$code), identity), use = "pairwise")
corrplot(cor_log_R_tot_z, diag = F, method = "ellipse", order = "hclust", 
         hclust.method = "average", addrect = 4, col = c1, tl.col = "black") 
mtext("Recruitment process error residuals", side = 3, line = 2.5, cex = 1.5)
# dev.off()
rm(cor_log_R_tot_z)

# Spawner abundance residuals
dev.new(width = 10, height = 10, mar = c(1,1,1,1))
# png(filename="log_S_tot_err_corr.png", width=10*0.9, height=10*0.9, units="in", res=200, type="cairo-png")
log_S_tot <- log(extract1(IPM_pp, "S_tot"))
log_S_tot_err <- colMeans(sweep(log_S_tot, 2, log(fish_data$S_tot_obs), "-"))
log_S_tot_err[!is.finite(log_S_tot_err)] <- NA
cor_log_S_tot_err <- cor(tapply(log_S_tot_err, list(fish_data$year, fish_data$code), identity), use = "pairwise")
corrplot(cor_log_S_tot_err, diag = F, method = "ellipse", order = "hclust", 
         hclust.method = "average", addrect = 3, col = c1, tl.col = "black") 
mtext("Spawner observation error residuals", side = 3, line = 2.5, cex = 1.5)
# dev.off()
rm(log_S_tot);rm(log_S_tot_err);rm(cor_log_S_tot_err);rm(c1)


#--------------------------------------------------------------------------------
# Time series of observed and fitted SPAWNER AGE DISTRIBUTIONS for each pop
# under hierarchical IPM
#--------------------------------------------------------------------------------

dev.new(width=16,height=10)
png(filename="q_fit_IPM.png", width=16*0.9, height=10*0.9, units="in", res=200, type="cairo-png")
par(mfrow=c(5,6), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))

q_IPM <- extract1(IPM_pp,"q")
q_obs <- fish_data[,grep("n_age", names(fish_data))]
q_obs <- sweep(q_obs, 1, ifelse(rowSums(q_obs) > 0, rowSums(q_obs), NA), "/")
init_NA <- fish_data$year
for(i in levels(fish_data$code))
  init_NA[fish_data$code==i] <- init_NA[fish_data$code==i] - min(init_NA[fish_data$code==i]) + 1
init_NA <- is.na(fish_data$S_tot_obs) & init_NA < 5

c1 <- "salmon"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.4)
c2 <- "dodgerblue4"
c2t <- col2rgb(c2)
c2t <- rgb(c2t[1], c2t[2], c2t[3], maxColorValue = 255, alpha = 255*0.4)

for(i in levels(fish_data$code))
{
  y1 <- fish_data$year[fish_data$code==i]
  plot(y1, q_obs[fish_data$code==i,1], pch = "",
       xlim = range(fish_data$year), ylim = c(0,1), 
       cex.axis = 1.2, las = 1, xlab = "", ylab = "")
  mtext(i, side = 3, line = 0.5, cex = par("cex")*1.5)
  if(i %in% levels(fish_data$code)[seq(1,29,6)]) mtext("P(spawner age)", side = 2, line = 3.5, cex = par("cex")*1.5)
  if(i %in% levels(fish_data$code)[25:29]) mtext("Year", side = 1, line = 3, cex = par("cex")*1.5)
  lines(y1, apply(q_IPM[,fish_data$code==i,1], 2, median), col = c1, lwd = 1.5)
  polygon(c(y1, rev(y1)), 
          c(apply(q_IPM[,fish_data$code==i,1], 2, quantile, 0.025), 
            rev(apply(q_IPM[,fish_data$code==i,1], 2, quantile, 0.975))),
          col = c1t, border = NA)
  lines(y1, apply(q_IPM[,fish_data$code==i,2], 2, median), col = c2, lwd = 1.5)
  polygon(c(y1, rev(y1)), 
          c(apply(q_IPM[,fish_data$code==i,2], 2, quantile, 0.025), 
            rev(apply(q_IPM[,fish_data$code==i,2], 2, quantile, 0.975))),
          col = c2t, border = NA)
  points(y1, q_obs[fish_data$code==i,1], pch=16, col = c1, cex = 1)
  points(y1, q_obs[fish_data$code==i,2], pch=16, col = c2, cex = 1)
}
plot(0:1, 0:1, pch = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
legend("topleft", c("3","4"), pch = 16, lwd = 1.5, col = c(c1,c2), cex = 1.5,
       title = "Spawner age (years)")

rm(list = c("q_IPM","q_obs","c1","c1t","c2","c2t","y1","init_NA"))
dev.off()


#--------------------------------------------------------------------------------
# Time series of observed and fitted MEAN SPAWNER AGE for each pop
# under hierarchical IPM
#--------------------------------------------------------------------------------

dev.new(width=16,height=10)
png(filename="mean_q_fit_IPM.png", width=16*0.9, height=10*0.9, units="in", res=200, type="cairo-png")
par(mfrow=c(5,6), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))

mean_q_IPM <- apply(extract1(IPM_pp,"q"), c(1,2), function(x) sum(x*3:5))
q_obs <- fish_data[,grep("n_age", names(fish_data))]
q_obs <- sweep(q_obs, 1, ifelse(rowSums(q_obs) > 0, rowSums(q_obs), NA), "/")
mean_q_obs <- rowSums(sweep(q_obs, 2, 3:5, "*"))
init_NA <- fish_data$year
for(i in levels(fish_data$code))
  init_NA[fish_data$code==i] <- init_NA[fish_data$code==i] - min(init_NA[fish_data$code==i]) + 1
init_NA <- is.na(fish_data$S_tot_obs) & init_NA < 5

c1 <- "blue4"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.4)

for(i in levels(fish_data$code))
{
  y1 <- fish_data$year[fish_data$code==i]
  plot(y1, mean_q_obs[fish_data$code==i], pch = "", 
       xlim = range(fish_data$year), 
       ylim = range(fish_data$mean_q_obs,
                    apply(mean_q_IPM[,!init_NA], 2, quantile, c(0.025,0.975)), na.rm = T), 
       cex.axis = 1.2, las = 1, xlab = "", ylab = "")
  mtext(i, side = 3, line = 0.5, cex = par("cex")*1.5)
  if(i %in% levels(fish_data$code)[seq(1,29,6)]) mtext("Mean spawner age", side = 2, line = 3.5, cex = par("cex")*1.5)
  if(i %in% levels(fish_data$code)[25:29]) mtext("Year", side = 1, line = 3, cex = par("cex")*1.5)
  lines(y1, apply(mean_q_IPM[,fish_data$code==i], 2, median), col = c1, lwd = 1.5)
  polygon(c(y1, rev(y1)), 
          c(apply(mean_q_IPM[,fish_data$code==i], 2, quantile, 0.025), 
            rev(apply(mean_q_IPM[,fish_data$code==i], 2, quantile, 0.975))),
          col = c1t, border = NA)
  points(y1, mean_q_obs[fish_data$code==i], pch=16, col = "black", cex = 1)
}

rm(list = c("mean_q_IPM","q_obs","mean_q_obs","c1","c1t","y1","init_NA"))
dev.off()


#--------------------------------------------------------------------------------
# Time series of fitted MEAN RECRUIT AGE for each pop under hierarchical IPM
# (data not shown b/c recruits not directly observed)
#--------------------------------------------------------------------------------

dev.new(width=16,height=10)
png(filename="mean_p_fit_IPM.png", width=16*0.9, height=10*0.9, units="in", res=200, type="cairo-png")
par(mfrow=c(5,6), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))

mean_p_IPM <- apply(extract1(IPM_pp,"p"), c(1,2), function(x) sum(x*3:5))

c1 <- "blue4"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.4)

for(i in levels(fish_data$code))
{
  y1 <- fish_data$year[fish_data$code==i]
  plot(y1, colMeans(mean_p_IPM[,fish_data$code==i]), pch = "", 
       xlim = range(fish_data$year), 
       ylim = range(apply(mean_p_IPM, 2, quantile, c(0.025,0.975))), 
       yaxs = "i", cex.axis = 1.2, las = 1, xlab = "", ylab = "")
  mtext(i, side = 3, line = 0.5, cex = par("cex")*1.5)
  if(i %in% levels(fish_data$code)[seq(1,29,6)]) mtext("Mean recruit age", side = 2, line = 3.5, cex = par("cex")*1.5)
  if(i %in% levels(fish_data$code)[25:29]) mtext("Year", side = 1, line = 3, cex = par("cex")*1.5)
  lines(y1, apply(mean_p_IPM[,fish_data$code==i], 2, median), col = c1, lwd = 2)
  polygon(c(y1, rev(y1)), 
          c(apply(mean_p_IPM[,fish_data$code==i], 2, quantile, 0.025), 
            rev(apply(mean_p_IPM[,fish_data$code==i], 2, quantile, 0.975))),
          col = c1t, border = NA)
}

rm(list = c("mean_p_IPM","c1","c1t","y1"))
dev.off()


#-----------------------------------------------------------------------------------------
# Cross-correlation among pops in 
# (1) Mean recruit age
# (2) Mean spawner age
# Correlation matrices are grouped by agglomerative clustering (UPGMA),
# with the largest number of clusters such that each cluster contains >= 3 pops
#-----------------------------------------------------------------------------------------

c1 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                             "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))(200)

# Mean recruit age
dev.new(width = 10, height = 10, mar = c(1,1,1,1))
# png(filename="mean_p_corr.png", width=10*0.9, height=10*0.9, units="in", res=200, type="cairo-png")
mean_p_IPM <- colMeans(apply(extract1(IPM_pp,"p"), c(1,2), function(x) sum(x*3:5)))
cor_mean_p <- cor(tapply(mean_p_IPM, list(fish_data$year, fish_data$code), identity), use = "pairwise")
corrplot(cor_mean_p, diag = F, method = "ellipse", order = "hclust", 
         hclust.method = "average", addrect = NULL, col = c1, tl.col = "black") 
mtext("Mean recruit age", side = 3, line = 2.5, cex = 1.5)
# dev.off()
rm(mean_p_IPM);rm(cor_mean_p)

dev.new(width = 10, height = 10, mar = c(1,1,1,1))
# png(filename="mean_q_corr.png", width=10*0.9, height=10*0.9, units="in", res=200, type="cairo-png")
mean_q_IPM <- colMeans(apply(extract1(IPM_pp,"q"), c(1,2), function(x) sum(x*3:5)))
cor_mean_q <- cor(tapply(mean_q_IPM, list(fish_data$year, fish_data$code), identity), use = "pairwise")
corrplot(cor_mean_q, diag = F, method = "ellipse", order = "hclust", 
         hclust.method = "average", addrect = NULL, col = c1, tl.col = "black") 
mtext("Mean spawner age", side = 3, line = 2.5, cex = 1.5)
# dev.off()
rm(mean_q_IPM);rm(cor_mean_q);rm(c1)


#--------------------------------------------------------------------------
# Shared brood-year productivity anomalies and process error innovations
#--------------------------------------------------------------------------

dev.new(width = 10, height = 7)
# png(filename="phi_timeseries.png", width=10, height=7, units="in", res=200, type="cairo-png")
par(mar = c(5.1,5.1,4.1,2.1))
mod <- IPM_pp
log_phi <- log(extract1(mod, "phi"))
log_phi_z_innov <- c(NA, stan_mean(mod,"log_phi_z")[-1])*stan_mean(mod, "sigma_log_phi")
y1 <- sort(unique(fish_data$year))
c1 <- "blue4"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.3)
c2 <- "green3"
c2t <- col2rgb(c2)
c2t <- rgb(c2t[1], c2t[2], c2t[3], maxColorValue = 255, alpha = 255*0.6)
plot(y1, colMeans(log_phi), pch = "", las = 1, cex.axis = 1.5, cex.lab = 1.8, 
     xaxt = "n", #log = "y", yaxt = "n",
     ylim = range(apply(log_phi, 2, quantile, c(0.025, 0.975)), log_phi_z_innov, na.rm = T),
     xlab = "Brood year", ylab = "Productivity anomaly")
axis(side = 1, at = y1[y1 %% 10 == 0], cex.axis = 1.5)
rug(y1[y1 %% 10 != 0], ticksize = -0.01)
rug(y1[y1 %% 10 != 0 & y1 %% 5 == 0], ticksize = -0.02)
# at <- maglab(10^par("usr")[3:4], log = T)
# axis(2, at$labat, cex.axis=1.2, las=1, labels = at$labat, cex.axis = 1.5)
abline(h = 0, lty = 2)
abline(v = y1, col = "lightgray")
lines(y1, log_phi_z_innov, col = c2t, lwd = 3)
lines(y1, colMeans(log_phi), type = "l", col = c1, lwd = 3)
polygon(c(y1, rev(y1)),
        c(apply(log_phi, 2, quantile, 0.025), rev(apply(log_phi, 2, quantile, 0.975))),
        col = c1t, border = NA)
legend("topright", c("anomalies","innovations"), lty = 1, lwd = 3, col = c(c1,c2t), bg = "white")
rm(list = c("mod","log_phi","y1","c1","c1t","log_phi_z_innov","c2","c2t"))
# dev.off()











