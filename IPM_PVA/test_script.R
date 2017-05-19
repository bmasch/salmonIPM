options(device = windows)
library(lme4)
library(rethinking)

#-----------------
# PLOTS
#-----------------

# Load data
fish_data <- read.table(file.path("~", "SalmonIPM", "IPM_PVA", "fish_data.txt"), sep = "\t", header = T)

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


#-----------------
# FIT MODELS
#-----------------

# Fit hierarchical spawner-recruit model to run reconstruction
RR_fit <- salmonIPM(fish_data = fish_data, model = "RR", pool_pops = TRUE, chains = 3, iter = 1000, warmup = 500,
                    control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))
print(RR_fit, pars = c("phi","R_hat"), include = FALSE)

# Fit hierarchical IPM
IPM_fit <- salmonIPM(fish_data = fish_data, model = "IPM", 
                     pars = c("mu_log_a","sigma_log_a","a",
                       "mu_log_b","sigma_log_b","b","rho_log_ab",
                       "sigma_log_phi","rho_log_phi","phi",
<<<<<<< HEAD
                       "mu_p","sigma_alr_p","gamma_alr_p","tau_alr_p","alr_p_z",
                       "sigma_proc","sigma_obs",
                       "S_tot","S_W_tot","S_H_tot","R_tot","log_R_tot_z"),
=======
                       "mu_p","sigma_gamma","R_gamma","gamma",
                       "sigma_alr_p","R_alr_p","p","sigma_proc","sigma_obs",
                       "S_tot","S_W_tot","S_H_tot","R_tot"),
>>>>>>> modstan
                     chains = 3, iter = 2000, warmup = 1000, 
                     control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

print(IPM_fit, pars = c("mu_log_a","sigma_log_a",
                        "mu_log_b","sigma_log_b","rho_log_ab",
                        "sigma_log_phi","rho_log_phi",
                        "mu_p","sigma_gamma","R_gamma","gamma",
                        "sigma_alr_p","R_alr_p","sigma_proc","sigma_obs"))
launch_shinystan(IPM_fit)

# Fit non-hierarchical IPM
IPM_fit3 <- salmonIPM(fish_data = fish_data, model = "IPM", pool_pops = FALSE,
                     pars = c("a","b","sigma_proc","rho_proc","sigma_proc","sigma_obs",
                              "gamma","sigma_alr_p","R_alr_p","p",
                              "S_tot","S_W_tot","S_H_tot","R_tot"),
                     chains = 3, iter = 2000, warmup = 1000, 
                     control = list(adapt_delta = 0.95, max_treedepth = 12))

print(IPM_fit3, pars = c("a","b","sigma_proc","rho_proc","sigma_proc","sigma_obs",
                         "gamma","sigma_alr_p","R_alr_p"))
launch_shinystan(IPM_fit3)



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
log_a <- log(extract1(IPM_fit2, "a"))
log_b <- log(extract1(IPM_fit2, "b"))
log_Rmax <- log_a - log_b

dev.new(width = 7, height = 7)
par(mfrow = c(2,2))

qqnorm(colMeans(log_b), main = "log(b [ha/spawner])")
qqline(colMeans(log_b))
qqnorm(colMeans(log_b)/tapply(dat$A, dat$pop, mean), main = "log(b [spawners^-1])")
qqline(colMeans(log_b)/tapply(dat$A, dat$pop, mean))

qqnorm(colMeans(log_Rmax), main = "log(Rmax [spawners/ha])")
qqline(colMeans(log_Rmax))
qqnorm(colMeans(log_Rmax)*tapply(dat$A, dat$pop, mean), main = "log(Rmax [spawners])")
qqline(colMeans(log_Rmax)*tapply(dat$A, dat$pop, mean))


