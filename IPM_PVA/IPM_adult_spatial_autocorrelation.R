library(rstan)
library(loo)
library(shinystan)
library(gtools)
library(lattice)
library(parallel)
library(corrplot)
source("stan_mean.r")
source("extract1.r")
source("stan_init.r")
source("run_recon.r")
source("dgnorm.r")
options(device=windows)


#===========================================================================
# DATA
#===========================================================================

# from Salmon Population Summary database, supp data from M. Ford and T. Cooney, habitat area from Damon Holzer
srchin_all <- read.table("srchin.4ss.all.txt", sep="\t", header=T)
srchin <- srchin_all
### TEMP: Impute one NA value of nS in Chamberlain 1986
srchin$nS[srchin$pop == "Chamberlain" & srchin$brood.yr == 1986] <- 
  mean(srchin$nS[srchin$pop == "Chamberlain"][1:5], na.rm = T)

# Indices of ocean conditions (monthly PDO and CUI)
# Pull out variables to try in the integrated model, based on an initial
# all-subsets selection fitting exercise using estimates of log(phi)
ocean <- read.table("ocean_conditions.txt", sep="\t", header=T)
ocean[,-1] <- scale(ocean[,-1]) # convert to std anomalies
ocean <- cbind(brood.yr = ocean$year - 2, ocean)
X <- ocean[match(sort(unique(srchin$brood.yr)), ocean$brood.yr),c("pdo.may","pdo.jun","pdo.jul")]
# X <- matrix(rowMeans(X), nrow = nrow(X))
X <- matrix(0, nrow(X))


# Data for STAN
stan_dat <- list(N = nrow(srchin),
                 pop = as.numeric(srchin$pop), year = as.numeric(factor(srchin$brood.yr)),
                 N_X = ncol(X), X = X,
                 N_pop_H = length(unique(srchin$pop[srchin$period %in% c("nonlocal","local")])),
                 which_pop_H = array(unique(as.numeric(srchin$pop)[srchin$period %in% c("nonlocal","local")]),
                                     dim = length(unique(as.numeric(srchin$pop)[srchin$period %in% c("nonlocal","local")]))),
                 N_S_obs = sum(!is.na(srchin$nS)),
                 which_S_obs = array(which(!is.na(srchin$nS)), dim = sum(!is.na(srchin$nS))),
                 S_tot_obs = replace(srchin$nS, is.na(srchin$nS) | srchin$nS==0, 1),
                 N_age = 3, max_age = 5,
                 n_age_obs = as.matrix(srchin[,c("n3","n4","n5")]),
                 N_H = sum(srchin$period %in% c("nonlocal","local")),
                 which_H = array(which(srchin$period %in% c("nonlocal","local")), 
                                 dim = max(sum(srchin$period %in% c("nonlocal","local")), 1)),
                 n_W_obs = array(srchin$nW[srchin$period %in% c("nonlocal","local")], 
                                 dim = max(sum(srchin$period %in% c("nonlocal","local")), 1)),
                 n_H_obs = array(srchin$nH[srchin$period %in% c("nonlocal","local")], 
                                 dim = max(sum(srchin$period %in% c("nonlocal","local")), 1)),
                 A = srchin$ha,
                 F_rate = srchin$hrate.w,
                 N_B = sum(srchin$wild.broodstk > 0),
                 which_B = array(which(srchin$wild.broodstk > 0),
                                 dim = max(sum(srchin$wild.broodstk > 0), 1)),
                 B_take_obs = srchin$wild.broodstk[srchin$wild.broodstk > 0])
if(stan_dat$N_pop_H == 0) stan_dat$which_pop_H <- array(1, dim = 1)
if(stan_dat$N_H == 0) 
{
  stan_dat$which_H <- array(1, dim = 1)
  stan_dat$n_W_obs <- array(1, dim = 1)
  stan_dat$n_H_obs <- array(1, dim = 1)
}

# Hydrologic distance matrix from Aimee Fullerton
hyd_dist <- read.csv("HydDist.csv", header = T, row.names = 1)
hyd_dist <- hyd_dist[dimnames(hyd_dist)[[1]] %in% levels(srchin$code), 
                     dimnames(hyd_dist)[[2]] %in% levels(srchin$code)]


#===========================================================================
# CALL STAN TO FIT MODELS
#===========================================================================

# Fit model!
stan_BH <- stan(file = "IPM_adult.stan",
                data = stan_dat, 
                init = stan_init(stan_dat, 3, fixedpop = FALSE), 
                pars = c("mu_log_a","sigma_log_a","a",
                         "mu_log_b","sigma_log_b","b",
                         "beta_log_phi","sigma_log_phi","rho_log_phi","phi",
                         "mu_p","sigma_alr_p","gamma_alr_p",
                         "mu_tau_alr_p","sigma_log_tau_alr_p","tau_alr_p","p",
                         "pHOS","B_rate_all",
                         "mu_sigma_proc","sigma_log_sigma_proc","sigma_proc","sigma_obs",
                         "S_tot","log_R_tot_z","R_tot","q"),
                chains = 3, iter = 2000, warmup = 1000, thin = 1, cores = 3,
                control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))


#===========================================================================
# PLOTS
#===========================================================================

# Recruitment process error residuals, log(R_tot) - log(R_tot_hat): cross-correlation b/w pops
dev.new(width = 10, height = 10, mar = c(1,1,1,1))
cor_log_R_tot_z <- cor(tapply(stan_mean(stan_BH,"log_R_tot_z"), list(stan_dat$year, stan_dat$pop), identity), use = "pairwise")
dimnames(cor_log_R_tot_z) <- list(srchin$code[match(levels(srchin$pop), srchin$pop)], 
                                  srchin$code[match(levels(srchin$pop), srchin$pop)])
corrplot(cor_log_R_tot_z, diag = F, method = "ellipse", order = "hclust", hclust.method = "average", 
         addrect = 5, tl.col = "black")


# Spatial correlogram of recruitment process error residuals
dev.new(width = 7, height = 7)
cor_log_R_tot_z <- cor(tapply(stan_mean(stan_BH,"log_R_tot_z"), list(stan_dat$year, stan_dat$pop), identity), use = "pairwise")
dimnames(cor_log_R_tot_z) <- list(srchin$code[match(levels(srchin$pop), srchin$pop)], 
                                  srchin$code[match(levels(srchin$pop), srchin$pop)])
cor_log_R_tot_z <- cor_log_R_tot_z[rownames(hyd_dist), colnames(hyd_dist)]
plot(hyd_dist[lower.tri(hyd_dist)], cor_log_R_tot_z[lower.tri(cor_log_R_tot_z)],
     xlab = "Hydrologic distance (km)", ylab = "Process error correlation",
     cex.axis = 1.2, cex.lab = 1.5, las = 1)
abline(h = 0, lty = 2)
  

