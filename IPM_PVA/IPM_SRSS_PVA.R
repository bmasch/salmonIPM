library(rstan)
library(loo)
library(shinystan)
library(gtools)
library(lattice)
library(parallel)
library(magicaxis)
source("stan_mean.r")
source("extract1.r")
source("IPM_adult_forecast.r")
source("step_ahead.r")
source("stan_init.r")
source("run_recon.r")
source("dgnorm.r")
source("fix_ff_filenames.r")
options(device=windows)


#===========================================================================
# DATA
#===========================================================================

srchin_all <- read.table("srchin.4ss.all.txt", sep="\t", header=T)
srchin <- srchin_all
### TEMP: Impute one NA value of nS in Chamberlain 1986
srchin$nS[srchin$pop == "Chamberlain" & srchin$brood.yr == 1986] <- 
  mean(srchin$nS[srchin$pop == "Chamberlain"][1:5], na.rm = T)

# Pad data with years through max_year
max_year <- max(srchin$brood.yr) + 50
year_aug <- sapply(tapply(srchin$brood.yr, srchin$pop, max), function(x) (x + 1):max_year)
pop_aug <- rep(names(year_aug), sapply(year_aug, length))
ha_aug <- rep(tapply(srchin$ha, srchin$pop, mean), times = sapply(year_aug, length))
srchin_aug <- data.frame(type = "future", pop = pop_aug, brood.yr = unlist(year_aug), 
                         ha = ha_aug, period = NA, nS = NA, n3 = 0, n4 = 0, n5 = 0,
                         p3 = NA, p4 = NA, p5 = NA, nH = 0, nW = 0, pHOS = 0, 
                         wild.broodstk = 0, hrate.w = 0,
                         row.names = NULL)
srchin <- rbind(data.frame(type = "past",
                           srchin[,c("pop","brood.yr","ha","period","nS","n3","n4","n5",
                                     "p3","p4","p5","nH","nW","pHOS","wild.broodstk","hrate.w")]), 
                srchin_aug)
srchin <- srchin[order(srchin$pop, srchin$brood.yr),]
row.names(srchin) <- NULL

X <- matrix(0, length(unique(srchin$brood.yr)))


#===========================================================================
# INTEGRATED MODEL
#===========================================================================

# Data for Stan
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

# Fit model!
stan_PVA <- stan(file = "IPM_adult.stan",
                data = stan_dat, 
                init = stan_init(stan_dat, chains = 3, fixedpop = FALSE), 
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


# # Write out results (total spawners and recruits per spawner by population)
# S_tot <- extract1(stan_PVA,"S_tot")[,srchin$brood.yr > max(srchin_all$brood.yr)]
# R_tot <- extract1(stan_PVA,"R_tot")[,srchin$brood.yr > max(srchin_all$brood.yr)]
# RS <- R_tot/S_tot
# pop <- srchin$pop[srchin$brood.yr > max(srchin_all$brood.yr)]
# brood_year <- srchin$brood.yr[srchin$brood.yr > max(srchin_all$brood.yr)]
# S_tot <- data.frame(pop = pop, brood_year = brood_year, t(S_tot))
# RS_tot <- data.frame(pop = pop, brood_year = brood_year, t(RS))
# write.table(S_tot, "S_tot.txt", sep = "\t", row.names = F)
# write.table(RS, "RS.txt", sep = "\t", row.names = F)
# rm(list = c("S_tot","R_tot","RS","pop","brood_year"))


#===========================================================================
# TRADITIONAL MULTILEVEL SPAWNER-RECRUIT MODEL
#===========================================================================

# Run reconstruction
srchin_rr <- run_recon(srchin)
rr_NA <- is.na(srchin$rr$R_tot) | is.na(srchin$nS)

# Data for Stan
stan_rr_dat <- list(N = sum(!reg_NA),
                     pop = as.numeric(srchin$pop)[!reg_NA], 
                     year = as.numeric(factor(srchin$brood.yr))[!reg_NA],
                     N_pop_H = length(unique(srchin$pop[srchin$period %in% c("nonlocal","local")])),
                     which_pop_H = array(unique(as.numeric(srchin$pop)[srchin$period %in% c("nonlocal","local")]),
                                         dim = length(unique(as.numeric(srchin$pop)[srchin$period %in% c("nonlocal","local")]))),
                     a_vary = 0, Rmax_vary = 0,
                     S_W = (replace(srchin$nS, srchin$nS == 0, 1)*(1 - srchin$pHOS))[!reg_NA],
                     S_H = (replace(srchin$nS, srchin$nS == 0, 1)*srchin$pHOS)[!reg_NA],
                     R = srchin_recon$R_tot[!reg_NA],
                     A = srchin$ha[!reg_NA])


# Function to generate initial values
stan_reg_init <- function() 
{
  if(exists("stan_dat")) 
    for(i in names(stan_reg_dat)) assign(i, stan_reg_dat[[i]])
  
  return(list(mu_log_a_W = runif(1,3,6), mu_log_a_H = runif(1,3,6), 
              sigma_log_a = runif(1,0.1,0.5),
              log_a_W_z = array(runif(max(pop),-1,1), dim = max(pop)), 
              log_a_H_z = array(runif(max(N_pop_H,1),-1,1), dim = max(N_pop_H,1)),
              mu_log_Rmax_W = runif(1,4,6), mu_log_Rmax_H = runif(1,4,6), 
              sigma_log_Rmax = runif(1,0.1,0.5),
              log_Rmax_W_z = array(runif(max(pop),-1,1), dim = max(pop)), 
              log_Rmax_H_z = array(runif(max(N_pop_H,1),-1,1), dim = max(N_pop_H,1)),
              sigma_log_phi = runif(1,0.1,0.5), 
              log_phi_z = array(rnorm(max(year),0,0.1), dim = max(year)),
              sigma = runif(1,0.1,2)))
}

# Fit model!
stan_LG_reg <- stan(file = "HW_STAN.stan",
                    data = c(stan_reg_dat), 
                    init = stan_reg_init, 
                    pars = c("mu_log_a_W","mu_log_a_H","sigma_log_a","a_W", "a_H", #"log_alpha",
                             "mu_log_Rmax_W","mu_log_Rmax_H","sigma_log_Rmax","Rmax_W","Rmax_H", #"log_rho",
                             "sigma_log_phi","phi","sigma","R_hat"),
                    chains = 3, iter = 10000, warmup = 2000, thin = 8, cores = 3,
                    control = list(adapt_delta = 0.8))



#===========================================================================
# FIGURES
#===========================================================================

# Time series of observed and fitted or predicted total spawners for each pop
dev.new(width=16,height=10)
# png(filename="IPM_SRSS_PVA_S_tot.png", width=16, height=10, units="in", res=200, type="cairo-png")
par(mfrow=c(4,6), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))
mod <- stan_PVA
S_tot <- extract1(mod,"S_tot")
c1 <- "darkgray"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.6)
c2 <- "blue"
c2t <- col2rgb(c2)
c2t <- rgb(c2t[1], c2t[2], c2t[3], maxColorValue = 255, alpha = 255*0.4)
iters <- 1:5
for(i in levels(srchin$pop))
{
  S_tot_i <- S_tot[,srchin$pop==i]
  S_tot_fit_i <- S_tot_i[,srchin$type[srchin$pop==i]=="past"]
  # S_tot_fit_pts_i <- density(log(S_tot_fit_i[,1]), 
  #                          from = min(apply(log(S_tot_i), 2, quantile, 0.025)), 
  #                          to = max(apply(log(S_tot_i), 2, quantile, 0.975)))$x
  # S_tot_fit_dens_i <- apply(log(S_tot_fit_i), 2, function(x) 
  #   density(x, from = min(apply(log(S_tot_i), 2, quantile, 0.025)), 
  #           to = max(apply(log(S_tot_i), 2, quantile, 0.975)))$y)
  S_tot_fore_i <- cbind(S_tot_fit_i[,ncol(S_tot_fit_i)],
                        S_tot_i[,srchin$type[srchin$pop==i]=="future"])
  # S_tot_fore_pts_i <- density(log(S_tot_fore_i[,1]), 
  #                            from = min(apply(log(S_tot_i), 2, quantile, 0.025)), 
  #                            to = max(apply(log(S_tot_i), 2, quantile, 0.975)))$x
  # S_tot_fore_dens_i <- apply(log(S_tot_fore_i), 2, function(x) 
  #   density(x, from = min(apply(log(S_tot_i), 2, quantile, 0.025)), 
  #           to = max(apply(log(S_tot_i), 2, quantile, 0.975)))$y)
  year_fit_i <- srchin$brood.yr[srchin$pop==i & srchin$type=="past"]
  year_fore_i <- c(year_fit_i[length(year_fit_i)], 
                   srchin$brood.yr[srchin$pop==i & srchin$type=="future"])
  
  plot(srchin$brood.yr[srchin$pop==i], srchin$nS[srchin$pop==i], pch="", cex=1.2, cex.axis=1.2, las=1,
       # ylim = c(0,max(apply(S_tot_i, 2, quantile, 0.975))), xlab="", ylab="")
       # ylim = c(min(apply(S_tot_i, 2, quantile, 0.025)),max(apply(S_tot_i, 2, quantile, 0.975))),
       ylim = range(c(range(S_tot_i[iters,]),
                      c(min(apply(S_tot_fit_i, 2, quantile, 0.025)),max(apply(S_tot_fit_i, 2, quantile, 0.975))),
                      range(replace(srchin$nS[srchin$pop==i], srchin$nS[srchin$pop==i]==0, 1), 
                            na.rm = T))),
       xlab="", ylab="", log = "y", yaxt = "n")
  at <- maglab(10^par("usr")[3:4], log = T)
  axis(2, at$labat, cex.axis=1.2, las=1,
       labels = sapply(log10(at$labat), function(i) as.expression(bquote(10^ .(i)))))
  # densregion(year_fit_i, exp(S_tot_fit_pts_i), t(S_tot_fit_dens_i), colmax = c1)
  lines(year_fit_i,apply(S_tot_fit_i, 2, quantile, 0.5), col=c1, lwd=3)
  polygon(c(year_fit_i, rev(year_fit_i)),
          c(apply(S_tot_fit_i, 2, quantile, 0.025), rev(apply(S_tot_fit_i, 2, quantile, 0.975))),
          col = c1t, border = NA)
  # densregion(year_fore_i, exp(S_tot_fore_pts_i), t(S_tot_fore_dens_i), colmax = c2)
  # lines(year_fore_i,apply(S_tot_fore_i, 2, quantile, 0.5), col=c2, lwd=3)
  # polygon(c(year_fore_i, rev(year_fore_i)),
  #         c(apply(S_tot_fore_i, 2, quantile, 0.025), rev(apply(S_tot_fore_i, 2, quantile, 0.975))),
  #         col = c2t, border = NA)
  for(j in iters)
  {
    lines(year_fit_i, S_tot_fit_i[j,], col = c1)
    lines(year_fore_i, S_tot_fore_i[j,], col = c2)
  }
  points(srchin$brood.yr[srchin$pop==i], srchin$nS[srchin$pop==i], type="b", pch=16)
  mtext(i, side=3, line=0.5, cex=1.2)
}
mtext("Year", outer = T, side=1, line=2, cex=1.2)
mtext("Spawners", outer = T, side=2, line=1.1, cex=1.2)
rm(list=c("mod","S_tot","S_tot_i","S_tot_fit_i","S_tot_fore_i","year_fit_i","year_fore_i",
          "c1","c1t","c2","c2t","at","iters"))
# dev.off()


# Time series of observed and fitted or predicted total spawners for 6 cherry-picked pops
dev.new(width=11,height=7)
# png(filename="IPM_SRSS_PVA_S_tot6.png", width=11, height=7, units="in", res=200, type="cairo-png")
par(mfrow=c(2,3), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))
pops <- c("BearValley","Camas","GrandRUpperMain","Lostine","MarshCR","Yankee")
mod <- stan_PVA
S_tot <- extract1(mod,"S_tot")
c1 <- "darkgray"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.6)
c2 <- "blue"
c2t <- col2rgb(c2)
c2t <- rgb(c2t[1], c2t[2], c2t[3], maxColorValue = 255, alpha = 255*0.4)
iters <- 1:5
for(i in pops)
{
  S_tot_i <- S_tot[,srchin$pop==i]
  S_tot_fit_i <- S_tot_i[,srchin$type[srchin$pop==i]=="past"]
  S_tot_fore_i <- cbind(S_tot_fit_i[,ncol(S_tot_fit_i)],
                        S_tot_i[,srchin$type[srchin$pop==i]=="future"])
  year_fit_i <- srchin$brood.yr[srchin$pop==i & srchin$type=="past"]
  year_fore_i <- c(year_fit_i[length(year_fit_i)], 
                   srchin$brood.yr[srchin$pop==i & srchin$type=="future"])
  
  plot(srchin$brood.yr[srchin$pop==i], srchin$nS[srchin$pop==i], pch="", cex=1.2, cex.axis=1.2, las=1,
       # ylim=c(min(apply(S_tot_i, 2, quantile, 0.025)),max(apply(S_tot_i, 2, quantile, 0.975))), 
       ylim = range(c(range(S_tot_i[iters,]),
                      c(min(apply(S_tot_fit_i, 2, quantile, 0.025)),max(apply(S_tot_fit_i, 2, quantile, 0.975))),
                      range(replace(srchin$nS[srchin$pop==i], srchin$nS[srchin$pop==i]==0, 1), 
                            na.rm = T))),
       xlab="", ylab="", log = "y", yaxt = "n")
  at <- maglab(10^par("usr")[3:4], log = T)
  axis(2, at$labat, cex.axis=1.2, las=1,
       labels = sapply(log10(at$labat), function(i) as.expression(bquote(10^ .(i)))))
  lines(year_fit_i,apply(S_tot_fit_i, 2, quantile, 0.5), col=c1, lwd=3)
  polygon(c(year_fit_i, rev(year_fit_i)),
          c(apply(S_tot_fit_i, 2, quantile, 0.025), rev(apply(S_tot_fit_i, 2, quantile, 0.975))),
          col = c1t, border = NA)
  # lines(year_fore_i,apply(S_tot_fore_i, 2, quantile, 0.5), col=c2, lwd=3)
  # polygon(c(year_fore_i, rev(year_fore_i)),
  #         c(apply(S_tot_fore_i, 2, quantile, 0.025), rev(apply(S_tot_fore_i, 2, quantile, 0.975))),
  #         col = c2t, border = NA)
  for(j in iters)
  {
    lines(year_fit_i, S_tot_fit_i[j,], col = c1)
    lines(year_fore_i, S_tot_fore_i[j,], col = c2)
  }
  points(srchin$brood.yr[srchin$pop==i], srchin$nS[srchin$pop==i], type="b", pch=16)
  mtext(i, side=3, line=0.5, cex=1.2)
}
mtext("Year", outer = T, side=1, line=2, cex=1.2)
mtext("Spawners", outer = T, side=2, line=1.1, cex=1.2)
rm(list=c("mod","pops","S_tot","S_tot_i","S_tot_fit_i","S_tot_fore_i","year_fit_i","year_fore_i",
          "c1","c1t","c2","c2t","at","iters"))
# dev.off()



# Time series of observed and fitted or predicted recruits per spawner for each pop
dev.new(width=16,height=10)
# png(filename="IPM_SRSS_PVA_RS.png", width=16, height=10, units="in", res=200, type="cairo-png")
par(mfrow=c(4,6), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))
mod <- stan_PVA
S_tot <- extract1(mod,"S_tot")
R_tot <- extract1(mod,"R_tot")
RS <- R_tot/S_tot
c1 <- "darkgray"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.6)
c2 <- "blue"
c2t <- col2rgb(c2)
c2t <- rgb(c2t[1], c2t[2], c2t[3], maxColorValue = 255, alpha = 255*0.4)
iters <- 1:5
for(i in levels(srchin$pop))
{
  RS_obs_i <- (srchin_recon$R_tot/srchin$nS)[srchin$pop==i]
  RS_i <- RS[,srchin$pop==i]
  RS_fit_i <- RS_i[,srchin$type[srchin$pop==i]=="past"]
  RS_fore_i <- cbind(RS_fit_i[,ncol(RS_fit_i)],
                        RS_i[,srchin$type[srchin$pop==i]=="future"])
  year_fit_i <- srchin$brood.yr[srchin$pop==i & srchin$type=="past"]
  year_fore_i <- c(year_fit_i[length(year_fit_i)], 
                   srchin$brood.yr[srchin$pop==i & srchin$type=="future"])
  
  plot(srchin$brood.yr[srchin$pop==i], RS_obs_i, pch = "", cex=1.2, cex.axis=1.2, las=1,
       # ylim=c(min(apply(RS_i, 2, quantile, 0.025)),max(apply(RS_i, 2, quantile, 0.975))), 
       ylim = range(c(range(RS_i[iters,]),
                      c(min(apply(RS_fit_i, 2, quantile, 0.025)),max(apply(RS_fit_i, 2, quantile, 0.975))),
                      range(replace(RS_obs_i, RS_obs_i==0 | is.infinite(RS_obs_i), NA), na.rm = T))),
       xlab="", ylab="", yaxt = "n", log = "y")
  at <- maglab(10^par("usr")[3:4], log = T)
  axis(2, at$labat, cex.axis=1.2, las=1,
       labels = at$labat)
  lines(year_fit_i,apply(RS_fit_i, 2, quantile, 0.5), col=c1, lwd=2)
  polygon(c(year_fit_i, rev(year_fit_i)),
          c(apply(RS_fit_i, 2, quantile, 0.025), rev(apply(RS_fit_i, 2, quantile, 0.975))),
          col = c1t, border = NA)
  # lines(year_fore_i,apply(RS_fore_i, 2, quantile, 0.5), col=c2, lwd=2)
  # polygon(c(year_fore_i, rev(year_fore_i)),
  #         c(apply(RS_fore_i, 2, quantile, 0.025), rev(apply(RS_fore_i, 2, quantile, 0.975))),
  #         col = c2t, border = NA)
  for(j in iters)
  {
    lines(year_fit_i, RS_fit_i[j,], col = c1)
    lines(year_fore_i, RS_fore_i[j,], col = c2)
  }
  points(srchin$brood.yr[srchin$pop==i], RS_obs_i, type = "b", pch = 16)
  mtext(i, side=3, line=0.5, cex=1.2)
}
mtext("Year", outer = T, side=1, line=2, cex=1.2)
mtext("Recruits / spawner", outer = T, side=2, line=1.1, cex=1.2)
rm(list=c("mod","S_tot","RS_obs_i","RS","RS_i","RS_fit_i","RS_fore_i",
          "year_fit_i","year_fore_i","c1","c1t","c2","c2t","at","iters"))
# dev.off()


# Time series of observed and fitted or predicted recruits per spawner for 6 cherry-picked pops
dev.new(width=11,height=7)
# png(filename="IPM_SRSS_PVA_RS6.png", width=11, height=7, units="in", res=200, type="cairo-png")
pops <- c("BearValley","Camas","GrandRUpperMain","Lostine","MarshCR","Yankee")
par(mfrow=c(2,3), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))
mod <- stan_PVA
S_tot <- extract1(mod,"S_tot")
R_tot <- extract1(mod,"R_tot")
RS <- R_tot/S_tot
c1 <- "darkgray"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.6)
c2 <- "blue"
c2t <- col2rgb(c2)
c2t <- rgb(c2t[1], c2t[2], c2t[3], maxColorValue = 255, alpha = 255*0.4)
iters <- 1:5
for(i in pops)
{
  RS_obs_i <- (srchin_recon$R_tot/srchin$nS)[srchin$pop==i]
  RS_i <- RS[,srchin$pop==i]
  RS_fit_i <- RS_i[,srchin$type[srchin$pop==i]=="past"]
  RS_fore_i <- cbind(RS_fit_i[,ncol(RS_fit_i)],
                     RS_i[,srchin$type[srchin$pop==i]=="future"])
  year_fit_i <- srchin$brood.yr[srchin$pop==i & srchin$type=="past"]
  year_fore_i <- c(year_fit_i[length(year_fit_i)], 
                   srchin$brood.yr[srchin$pop==i & srchin$type=="future"])
  
  plot(srchin$brood.yr[srchin$pop==i], RS_obs_i, pch = "", cex=1.2, cex.axis=1.2, las=1,
       # ylim=c(min(apply(RS_fit_i, 2, quantile, 0.025)),max(apply(RS_fit_i, 2, quantile, 0.975))), 
       ylim = range(c(range(RS_i[iters,]),
                      c(min(apply(RS_fit_i, 2, quantile, 0.025)),max(apply(RS_fit_i, 2, quantile, 0.975))),
                      range(replace(RS_obs_i, RS_obs_i==0 | is.infinite(RS_obs_i), NA), na.rm = T))),
       xlab="", ylab="", yaxt = "n", log = "y")
  at <- maglab(10^par("usr")[3:4], log = T)
  axis(2, at$labat, cex.axis=1.2, las=1, labels = at$labat)
  lines(year_fit_i,apply(RS_fit_i, 2, quantile, 0.5), col=c1, lwd=2)
  polygon(c(year_fit_i, rev(year_fit_i)),
          c(apply(RS_fit_i, 2, quantile, 0.025), rev(apply(RS_fit_i, 2, quantile, 0.975))),
          col = c1t, border = NA)
  # lines(year_fore_i,apply(RS_fore_i, 2, quantile, 0.5), col=c2, lwd=2)
  # polygon(c(year_fore_i, rev(year_fore_i)),
  #         c(apply(RS_fore_i, 2, quantile, 0.025), rev(apply(RS_fore_i, 2, quantile, 0.975))),
  #         col = c2t, border = NA)
  for(j in iters)
  {
    lines(year_fit_i, RS_fit_i[j,], col = c1)
    lines(year_fore_i, RS_fore_i[j,], col = c2)
  }
  points(srchin$brood.yr[srchin$pop==i], RS_obs_i, type = "b", pch = 16)
  mtext(i, side=3, line=0.5, cex=1.2)
}
mtext("Year", outer = T, side=1, line=2, cex=1.2)
mtext("Recruits / spawner", outer = T, side=2, line=1.1, cex=1.2)
rm(list=c("mod","pops","S_tot","RS_obs_i","RS","RS_i","RS_fit_i","RS_fore_i",
          "year_fit_i","year_fore_i","c1","c1t","c2","c2t","at","iters"))
# dev.off()


# Shared brood-year productivity anomalies
dev.new(height = 10, width = 10)
# png(filename="IPM_SRSS_PVA_year_effects.png", width=10, height=10, units="in", res=200, type="cairo-png")
par(mar = c(5.1,5.1,4.1,2.1))
mod <- stan_PVA
phi <- extract1(mod, "phi")
phi_fit <- phi[,1:max(stan_dat$year[srchin$type=="past"])]
phi_fore <- phi[,(max(stan_dat$year[srchin$type=="past"]) + 1):ncol(phi)]
phi_fore <- cbind(phi_fit[,ncol(phi_fit)], phi_fore)
years <- sort(unique(srchin$brood.yr))
years_fit <- years[1:max(stan_dat$year[srchin$type=="past"])]
years_fore <- years[(max(stan_dat$year[srchin$type=="past"]) + 1):length(years)]
years_fore <- c(years_fit[length(years_fit)], years_fore)
iters <- 1:5
c1 <- "darkgray"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.6)
c2 <- "blue"
c2t <- col2rgb(c2)
c2t <- rgb(c2t[1], c2t[2], c2t[3], maxColorValue = 255, alpha = 255*0.4)
plot(years, colMeans(phi), pch = "", las = 1, cex.axis = 1.5, cex.lab = 1.8, log = "y", yaxt = "n",
     ylim = c(min(apply(phi, 2, quantile, 0.025)), max(apply(phi, 2, quantile, 0.975))),
     xlab = "Brood year", ylab = "Productivity anomaly")
at <- maglab(10^par("usr")[3:4], log = T)
axis(2, at$labat, cex.axis=1.2, las=1, labels = at$labat, cex.axis = 1.5)
lines(years_fit, colMeans(phi_fit), col = c1, lwd = 3)
polygon(c(years_fit, rev(years_fit)),
        c(apply(phi_fit, 2, quantile, 0.025), rev(apply(phi_fit, 2, quantile, 0.975))),
        col = c1t, border = NA)
# lines(years_fore, colMeans(phi_fore), col = c2, lwd = 3)
# polygon(c(years_fore, rev(years_fore)),
#         c(apply(phi_fore, 2, quantile, 0.025), rev(apply(phi_fore, 2, quantile, 0.975))),
#         col = c2t, border = NA)
for(i in iters)
{
  lines(years_fit, phi_fit[i,], col = c1t)
  lines(years_fore, phi_fore[i,], col = c2)
}
rm(list = c("mod","phi","phi_fit","phi_fore","years","years_fit","years_fore",
            "c1","c1t","c2","c2t","iters"))
# dev.off()



# Probability of quasi-extinction by population
dev.new(height = 10, width = 10)
# png(filename="IPM_SRSS_PVA_PQE.png", width=10, height=10, units="in", res=200, type="cairo-png")
par(oma = c(0,5,0,0))
qet <- 50     # set quasi-extinction threshold (4-yr moving average)
pop <- srchin$pop[srchin$type=="future"]
mod <- stan_PVA
S_tot <- t(extract1(mod, "S_tot")[,srchin$type=="future"])
pqe <- aggregate(S_tot, list(pop = pop), function(x) any(rollmean(x, 4) < qet))
pqe <- data.frame(pop = pqe[,1], pqe = rowMeans(pqe[,-1]))
pqe <- pqe[order(pqe$pqe),]

barplot(pqe$pqe, names.arg = pqe$pop, horiz = TRUE, 
        las = 1, cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5,
        xlab = "Probability of quasi-extinction", 
        main = paste("50-year quasi-extinction risk \n QET = ", qet, " spawners (4-year moving average)", sep = ""))
rm(list=c("qet","pop","mod","S_tot","pqe"))
# dev.off()


# CDFs of probability of quasi-extinction by population
dev.new(width=16,height=10)
# png(filename="IPM_SRSS_PVA_PQE_CDF.png", width=16, height=10, units="in", res=200, type="cairo-png")
par(mfrow=c(4,6), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))
pop <- srchin$pop[srchin$type=="future"]
mod <- stan_PVA
S_tot <- t(extract1(mod, "S_tot")[,srchin$type=="future"])

for(i in levels(srchin$pop))
{
  S_tot_ma <- apply(S_tot[pop==i,], 2, rollmean, 4)
  S_tot_ma_min <- apply(S_tot_ma, 2, min)
  plot.ecdf(S_tot_ma_min, xlab = "", ylab = "", main = i, xlim = c(0, max(S_tot_ma_min)),
            xaxs = "i", las = 1, cex.axis = 1.2, cex.lab = 1.5)
  segments(c(50, 50), c(0, ecdf(S_tot_ma_min)(50)), 
           c(50, 0), c(ecdf(S_tot_ma_min)(50), ecdf(S_tot_ma_min)(50)), col = "red")
}
mtext("QET", outer = T, side=1, line=2, cex=1.2)
mtext("Probability of quasi-extinction", outer = T, side=2, line=1.1, cex=1.2)
rm(list = c("pop","mod","S_tot","S_tot_ma","S_tot_ma_min"))
# dev.off()




