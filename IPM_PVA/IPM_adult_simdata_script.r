#library(rstan)
#library(loo)
#library(shinystan)
#library(gtools)
#library(lattice)
#library(parallel)
#library(zoo)
#library(magicaxis)
#source("stan_mean.r")
# source("extract1.r")
# source("IPM_adult_sim.r")
# source("stan_init.r")
# source("run_recon.r")
# source("dgnorm.r")
# options(device=windows)


#===========================================================================
# SIMULATE DATA
#===========================================================================

# Simulate data
pop <- rep(1:12, each = 50)
year <- rep(1:50, 12)
sim_dat <- IPM_adult_sim(pars = list(mu_log_a = 2, sigma_log_a = 0.5,
                                     mu_log_b = -1, sigma_log_b = 0.5,
                                     beta_log_phi = 1, rho_log_phi = 0.7, sigma_log_phi = 0.5, 
                                     mu_sigma_proc = 0.3, sigma_log_sigma_proc = 0.5, sigma_obs = 0.5, 
                                     mu_p = c(0.05, 0.55, 0.4), sigma_alr_p = c(0.5, 0.5), 
                                     mu_tau_alr_p = c(0.5, 0.5), sigma_log_tau_alr_p = c(0.3, 0.3)),
                         pop = pop, year = year, X = matrix(0, max(year), 1),
                         N_age = 3, max_age = 5, A = 100,
                         S_H_tot = replace(sample(1:1000, length(year), replace = T), year <= 30, 0),
                         F_rate = rbeta(length(pop), 7, 3), 
                         B_rate = replace(runif(length(pop), 0.01, 0.1), year <= 30, 0),
                         n_age_tot_obs = 50, n_HW_tot_obs = 50)

stan_dat <- list(N = sim_dat$sim_dat$N,
                 pop = pop, year = year,
                 N_X = ncol(sim_dat$sim_dat$X), X = sim_dat$sim_dat$X,
                 N_pop_H = length(unique(pop[sim_dat$pars_out$pHOS > 0])),
                 which_pop_H = array(unique(pop[sim_dat$pars_out$pHOS > 0]),
                                     dim = length(unique(pop[sim_dat$pars_out$pHOS > 0]))),
                 N_S_obs = sum(!is.na(sim_dat$sim_dat$S_tot_obs)),
                 which_S_obs = array(which(!is.na(sim_dat$sim_dat$S_tot_obs)), dim = sum(!is.na(sim_dat$sim_dat$S_tot_obs))),
                 S_tot_obs = replace(sim_dat$sim_dat$S_tot_obs, is.na(sim_dat$sim_dat$S_tot_obs), 1),
                 N_age = 3, max_age = 5,
                 n_age_obs = sim_dat$sim_dat$n_age_obs,
                 N_H = sum(sim_dat$pars_out$pHOS > 0),
                 which_H = array(which(sim_dat$pars_out$pHOS > 0), dim = sum(sim_dat$pars_out$pHOS > 0)),
                 n_W_obs = array(sim_dat$sim_dat$n_W_obs[sim_dat$pars_out$pHOS > 0],
                                 dim = sum(sim_dat$pars_out$pHOS > 0)),
                 n_H_obs = array(sim_dat$sim_dat$n_H_obs[sim_dat$pars_out$pHOS > 0],
                                 dim = sum(sim_dat$pars_out$pHOS > 0)),
                 A = rep(sim_dat$sim_dat$A, sim_dat$sim_dat$N),
                 F_rate = sim_dat$sim_dat$F_rate,
                 N_B = sum(sim_dat$sim_dat$B_take > 0),
                 which_B = array(which(sim_dat$sim_dat$B_take > 0), dim = sum(sim_dat$sim_dat$B_take > 0)),
                 B_take_obs = sim_dat$sim_dat$B_take[sim_dat$sim_dat$B_take > 0])
if(stan_dat$N_pop_H == 0) stan_dat$which_pop_H <- array(1, dim = 1)
if(stan_dat$N_H == 0)
{
  stan_dat$which_H <- array(1, dim = 1)
  stan_dat$n_W_obs <- array(1, dim = 1)
  stan_dat$n_H_obs <- array(1, dim = 1)
}


#===========================================================================
# CALL STAN TO FIT MODELS
#===========================================================================

# Partial pooling across populations
fixedpop <- FALSE
stan_BH <- stan(file = ifelse(fixedpop, "IPM_adult_fixedpop.stan", "IPM_adult.stan"),
                data = stan_dat, 
                init = stan_init(stan_dat, 3, fixedpop), 
                pars = c("mu_log_a","sigma_log_a","a",
                         "mu_log_b","sigma_log_b","b",
                         "beta_log_phi","sigma_log_phi","rho_log_phi","phi",
                         "mu_p","sigma_alr_p","gamma_alr_p",
                         "mu_tau_alr_p","sigma_log_tau_alr_p","tau_alr_p","p",
                         "pHOS","B_rate_all",
                         "mu_sigma_proc","sigma_log_sigma_proc","sigma_proc","sigma_obs",
                         "S_tot","R_tot","q"),
                chains = 3, iter = 2000, warmup = 1000, thin = 1, cores = 3,
                control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

# No pooling across populations
fixedpop <- TRUE
stan_BH_fixedpop <- stan(file = ifelse(fixedpop, "IPM_adult_fixedpop.stan", "IPM_adult.stan"),
                         data = stan_dat, 
                         init = stan_init(stan_dat, 3, fixedpop), 
                         pars = c("a","b","beta_log_proc","rho_proc","sigma_proc","gamma_p_arr",
                                  "tau_alr_p","p","pHOS","B_rate_all","sigma_obs","S_tot","R_tot","q"),
                         chains = 3, iter = 2000, warmup = 1000, thin = 1, cores = 3,
                         control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))


#===========================================================================
# PLOTS
#===========================================================================

# Time series of observed and estimated S_tot under unpooled and partially pooled models, 
# panel for each pop
# png(filename="S_tot_simdata.png", width=16*0.6, height=10*0.6, units="in", res=200, type="cairo-png")
dev.new(width=16,height=10)
par(mfrow=c(3,4), mar=c(1,2,2,1), oma=c(4.1,4.1,0,0))
S_tot <- extract1(stan_BH,"S_tot")
S_tot_fixedpop <- extract1(stan_BH_fixedpop,"S_tot")
c1 <- "darkgray"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.6)
for(i in unique(pop))
{
  plot(year[pop==i], sim_dat$sim_dat$S_tot_obs[pop==i], pch="", cex=1.2, cex.axis=1.2, las=1,
       ylim=c(min(apply(S_tot[,pop==i], 2, quantile, 0.025)),
              max(apply(S_tot[,pop==i], 2, quantile, 0.975))), 
       xlab="", ylab="", log = "y", yaxt = "n")
  at <- maglab(10^par("usr")[3:4], log = T)
  axis(2, at$labat, cex.axis=1.2, las=1,
       labels = sapply(log10(at$labat), function(i) as.expression(bquote(10^ .(i)))))
  lines(year[pop==i], apply(S_tot[,pop==i], 2, quantile, 0.5), lwd=2)
  lines(year[pop==i], apply(S_tot[,pop==i], 2, quantile, 0.025), lwd=1)
  lines(year[pop==i], apply(S_tot[,pop==i], 2, quantile, 0.975), lwd=1)
  lines(year[pop==i], apply(S_tot_fixedpop[,pop==i], 2, quantile, 0.5), col=c1, lwd=2)
  lines(year[pop==i], apply(S_tot_fixedpop[,pop==i], 2, quantile, 0.025), col=c1, lwd=1)
  lines(year[pop==i], apply(S_tot_fixedpop[,pop==i], 2, quantile, 0.975), col=c1, lwd=1)
  # polygon(c(year[pop==i], rev(year[pop==i])),
  #         c(apply(S_tot[,pop==i], 2, quantile, 0.025), rev(apply(S_tot[,pop==i], 2, quantile, 0.975))),
  #         col = c1t, border = NA)
  points(year[pop==i],  sim_dat$sim_dat$S_tot_obs[pop==i], pch = 16, type = "p")
}
mtext("Year", outer = T, side=1, line=2, cex=2*par("cex"))
mtext("Spawners", outer = T, side=2, line=1.1, cex=2*par("cex"))
rm(list=c("S_tot","S_tot_fixedpop","c1","c1t"))
# dev.off()


# Posterior densities under partially pooled model vs. true values
# png(filename="posterior_densities_simdata.png", width=16*0.6, height=10*0.6, units="in", res=200, type="cairo-png")
dev.new(width=16, height=10)
par(mfrow=c(3,4), mar=c(4.1,2,2,1), oma=c(0,4.1,0,0))
mod <- stan_BH
plot(density(extract1(mod,"mu_log_a")), xlab = bquote(mu[a]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$mu_log_a, lwd = 3)
plot(density(extract1(mod,"sigma_log_a")), xlab = bquote(sigma[a]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$sigma_log_a, lwd = 3)
plot(density(extract1(mod,"mu_log_b")), xlab = bquote(mu[b]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$mu_log_b, lwd = 3)
plot(density(extract1(mod,"sigma_log_b")), xlab = bquote(sigma[b]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$sigma_log_b, lwd = 3)
plot(density(extract1(mod,"beta_log_phi")), xlab = bquote(beta[phi]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$beta_log_phi, lwd = 3)
plot(density(extract1(mod,"rho_log_phi")), xlab = bquote(rho[phi]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = plogis(sim_dat$pars_out$logit_rho_log_phi), lwd = 3)
plot(density(extract1(mod,"sigma_log_phi")), xlab = bquote(sigma[phi]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$sigma_log_phi, lwd = 3)
plot(density(extract1(mod,"mu_p")[,1]), xlab = bquote(mu[p]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1, 
     xlim = range(extract1(mod,"mu_p")), col = "lightgray")
abline(v = sim_dat$pars_out$mu_p[1], lwd = 3, col = "lightgray")
dd <- density(extract1(mod,"mu_p")[,2])
lines(dd$x, dd$y, col = "darkgray")
abline(v = sim_dat$pars_out$mu_p[2], lwd = 3, col = "darkgray")
dd <- density(extract1(mod,"mu_p")[,3])
lines(dd$x, dd$y)
abline(v = sim_dat$pars_out$mu_p[3], lwd = 3)
plot(density(extract1(mod,"sigma_alr_p")[,1]), xlab = bquote(sigma[p]), ylab="", main="", 
     cex.axis=1.2, cex.lab=2, las=1, col = "darkgray")
abline(v = sim_dat$pars_out$sigma_alr_p[1], lwd = 3, col = "darkgray")
dd <- density(extract1(mod,"sigma_alr_p")[,2])
lines(dd$x, dd$y)
abline(v = sim_dat$pars_out$sigma_alr_p[2], lwd = 3)
plot(density(extract1(mod,"mu_sigma_proc")), xlab = bquote(mu[sigma[proc]]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$mu_sigma_proc, lwd = 3)
plot(density(extract1(mod,"sigma_log_sigma_proc")), xlab = bquote(sigma[sigma[proc]]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$sigma_log_sigma_proc, lwd = 3)
plot(density(extract1(mod,"sigma_obs")), xlab = bquote(sigma[obs]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$sigma_obs, lwd = 3)
mtext("Posterior density", side = 2, line = 2, outer = T, cex = 2*par("cex"))
rm(list=c("mod","dd"))
# dev.off()


# Posterior densities of pop-level random effects under unpooled and partially pooled models
# png(filename="random_effects_simdata.png", width=10*0.8, height=10*0.8, units="in", res=200, type="cairo-png")
dev.new(width=10, height=10)
par(mfcol=c(2,1), mar=c(2,5,1,1), oma=c(3.1,0,0,0))
log_a <- log(extract1(stan_BH,"a"))
log_a_fixedpop <- log(extract1(stan_BH_fixedpop,"a"))
log_b <- log(extract1(stan_BH,"b"))
log_b_fixedpop <- log(extract1(stan_BH_fixedpop,"b"))
plot(c(1,ncol(log_a)), range(apply(cbind(log_a, log_a_fixedpop), 2, quantile, c(0.01,0.99))), 
     pch="", ylab = bquote(log(italic(a))), xlab = "", cex.lab=1.5, cex.axis=1.2, las = 1)
for(i in 1:ncol(log_a))
{
  dd1 <- density(log_a_fixedpop[,i])
  # lines(i - 0.4*dd1$y/max(dd1$y), dd1$x, col = "darkgray")
  polygon(c(i - 0.4*dd1$y/max(dd1$y), i, i), c(dd1$x, range(dd1$x)),
          col = "darkgray", border = NA)
  dd2 <- density(log_a[,i])
  # lines(i + 0.4*dd2$y/max(dd2$y), dd2$x)
  polygon(c(i + 0.4*dd2$y/max(dd2$y), i, i), c(dd2$x, range(dd2$x)),
          col = "black", border = NA)
  points(i, log(sim_dat$pars_out$a[i]), pch = 16, col = "red")
}
legend("topleft", c("single-population","multi-population","true value"), pch=c(15,15,16),
       col=c("darkgray","black","red"))
plot(c(1,ncol(log_b)), range(apply(cbind(log_b, log_b_fixedpop), 2, quantile, c(0.01,0.99))), 
     pch="", ylab = bquote(log(italic(b))), xlab = "", cex.lab=1.5, cex.axis=1.2, las = 1)
for(i in 1:ncol(log_b))
{
  dd1 <- density(log_b_fixedpop[,i])
  # lines(i - 0.4*dd1$y/max(dd1$y), dd1$x, col = "darkgray")
  polygon(c(i - 0.4*dd1$y/max(dd1$y), i, i), c(dd1$x, range(dd1$x)),
          col = "darkgray", border = NA)
  dd2 <- density(log_b[,i])
  # lines(i + 0.4*dd2$y/max(dd2$y), dd2$x)
  polygon(c(i + 0.4*dd2$y/max(dd2$y), i, i), c(dd2$x, range(dd2$x)),
          col = "black", border = NA)
  points(i, log(sim_dat$pars_out$b[i]), pch = 16, col = "red")
}
mtext("Population", side=1, outer=T, line=1, cex=1.5*par("cex"))
rm(list=c("log_a","log_a_fixedpop","log_b","log_b_fixedpop","dd1","dd2"))
# dev.off()











