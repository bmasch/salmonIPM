#library(rstan)
#library(loo)
#library(shinystan)
#library(gtools)
#library(lattice)
#library(parallel)
#source("stan_mean.r")
#source("extract1.r")
#source("stan_init.r")
#source("run_recon.r")
#source("dgnorm.r")
#options(device=windows)


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
                         "S_tot","R_tot","q"),
                chains = 3, iter = 2000, warmup = 1000, thin = 1, cores = 3,
                control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

stan_BH_fixedpop <- stan(file = "IPM_adult_fixedpop.stan",
                         data = stan_dat, 
                         init = stan_init(stan_dat, 3, fixedpop = TRUE), 
                         pars = c("a","b","beta_proc","rho_proc","sigma_proc","gamma_p_arr","tau_alr_p","p",
                                  "p_HOS","B_rate_all","sigma_obs","S_tot","R_tot","q"),
                         # chains = 3, iter = 2000, warmup = 1000, thin = 1, cores = 3,
                         chains = 3, iter = 1500, warmup = 1000, thin = 1, cores = 3,
                         control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))




#===========================================================================
# TRADITIONAL SPAWNER-RECRUIT ANALYSIS BY MULTILEVEL REGRESSION
#===========================================================================

# Data for STAN
srchin_recon <- run_recon(srchin)
reg_NA <- is.na(srchin_recon$R_tot) | is.na(srchin$nS)

stan_reg_dat <- list(N = sum(!reg_NA),
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
  if(exists("stan_reg_dat")) 
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
stan_BH_reg <- stan(file = "HW_STAN.stan",
                    data = c(stan_reg_dat), 
                    init = stan_reg_init, 
                    pars = c("mu_log_a_W","mu_log_a_H","sigma_log_a","a_W", "a_H", #"log_alpha",
                             "mu_log_Rmax_W","mu_log_Rmax_H","sigma_log_Rmax","Rmax_W","Rmax_H", #"log_rho",
                             "sigma_log_phi","phi","sigma","R_hat"),
                    chains = 3, iter = 10000, warmup = 2000, thin = 8, cores = 3,
                    control = list(adapt_delta = 0.8))



#===========================================================================
# PLOTS
#===========================================================================

# Time series of observed and estimated S_tot, panel for each pop
dev.new(width=16,height=10)
par(mfrow=c(4,6), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))
mod <- stan_BH
S_tot <- extract(mod,"S_tot")[[1]]/1000
for(i in levels(srchin$pop))
{
  plot(srchin$brood.yr[srchin$pop==i], srchin$nS[srchin$pop==i]/1000, type="b", pch=16, cex=1.2, cex.axis=1.2, las=1,
       ylim=c(0,max(apply(S_tot[,srchin$pop==i], 2, quantile, 0.975))), xlab="", ylab="")
  lines(srchin$brood.yr[srchin$pop==i],apply(S_tot[,srchin$pop==i], 2, quantile, 0.025), col="blue")
  lines(srchin$brood.yr[srchin$pop==i],apply(S_tot[,srchin$pop==i], 2, quantile, 0.975), col="blue")
  lines(srchin$brood.yr[srchin$pop==i],apply(S_tot[,srchin$pop==i], 2, quantile, 0.5), col="blue", lwd=2)
  #   lines(srchin$brood.yr[srchin$pop==i],colMeans(S_tot[,srchin$pop==i]), col="blue", lwd=2)
  mtext(i, side=3, line=0.5, cex=1.2)
  if(i %in% levels(srchin$pop)[19:24]) mtext("Year", side=1, line=3, cex=1.2)
  if(i %in% levels(srchin$pop)[seq(1,24,6)]) mtext("Spawners (1000s)", side=2, line=3.1, cex=1.2)
}
rm(list=c("mod","S_tot"))


# Time series of observed and estimated S_tot, panel for each of 6 cherry-picked pops
# png(filename="S_tot.png", width=11, height=7, units="in", res=200, type="cairo-png")
dev.new(width=11,height=7)
par(mfrow=c(2,3), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))
mod <- stan_BH
S_tot <- extract(mod,"S_tot")[[1]]/1000
pops <- c("Lostine","Tucannon","MarshCR","SFSEast","Pahsimeroi","SREastFork")
for(i in pops)
{
  bb <- "blue4"
  plot(srchin$brood.yr[srchin$pop==i], srchin$nS[srchin$pop==i]/1000, type="b", pch=16, cex=1.2, cex.axis=1.2, las=1,
       ylim=c(0,max(apply(S_tot[,srchin$pop==i], 2, quantile, 0.975))), xlab="", ylab="")
  lines(srchin$brood.yr[srchin$pop==i],apply(S_tot[,srchin$pop==i], 2, quantile, 0.5), col=bb, lwd=3)
  bb <- col2rgb(bb)
  bb <- rgb(bb[1], bb[2], bb[3], maxColorValue=255, alpha=255*0.3)
  polygon(c(srchin$brood.yr[srchin$pop==i], rev(srchin$brood.yr[srchin$pop==i])), 
          c(apply(S_tot[,srchin$pop==i], 2, quantile, 0.025), 
            rev(apply(S_tot[,srchin$pop==i], 2, quantile, 0.955))),
          col = bb, border = NA)
  #   lines(srchin$brood.yr[srchin$pop==i],colMeans(S_tot[,srchin$pop==i]), col="blue4", lwd=3)
  mtext(i, side=3, line=0.5, cex=1.2)
  if(i %in% pops[4:6]) mtext("Year", side=1, line=3, cex=1.2)
  if(i %in% pops[c(1,4)]) mtext("Spawners (1000s)", side=2, line=3.1, cex=1.2)
}
rm(list=c("mod","pops","S_tot","bb"))
# dev.off()


# Time series of observed and estimated pHOS, panel for each pop
dev.new(width=15,height=9)
par(mfrow=c(3,5), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))
mod <- stan_BH
pHOS <- extract(mod,"pHOS")[[1]]
for(i in sort(levels(srchin$pop)[stan_dat$which_pop_H]))
{
  plot(srchin$brood.yr[srchin$pop==i], srchin$pHOS[srchin$pop==i], 
       type="b", pch=16, cex=1.2, cex.axis=1.2, las=1, ylim=c(0,1), xlab="", ylab="")
  lines(srchin$brood.yr[stan_dat$which_H][srchin$pop[stan_dat$which_H]==i],
        apply(pHOS[,srchin$pop[stan_dat$which_H]==i], 2, quantile, 0.025), col="blue")
  lines(srchin$brood.yr[stan_dat$which_H][srchin$pop[stan_dat$which_H]==i],
        apply(pHOS[,srchin$pop[stan_dat$which_H]==i], 2, quantile, 0.975), col="blue")
  lines(srchin$brood.yr[stan_dat$which_H][srchin$pop[stan_dat$which_H]==i],
        apply(pHOS[,srchin$pop[stan_dat$which_H]==i], 2, quantile, 0.5), col="blue", lwd=2)
  #   lines(srchin$brood.yr[srchin$pop==i],colMeans(S_tot[,srchin$pop==i]), col="blue", lwd=2)
  mtext(i, side=3, line=0.5, cex=1.2)
  if(i %in% sort(levels(srchin$pop)[stan_dat$which_pop_H])[11:14]) mtext("Year", side=1, line=3, cex=1.2)
  if(i %in% sort(levels(srchin$pop)[stan_dat$which_pop_H])[seq(1,14,5)]) mtext(bquote(p[HOS]), side=2, line=3.1, cex=1.2)
}
rm(list=c("mod","pHOS"))


# Time series of observed and estimated pHOS, panel for each of 6 cherry-picked pops
# png(filename="pHOS.png", width=11, height=7, units="in", res=200, type="cairo-png")
dev.new(width=11,height=7)
par(mfrow=c(2,3), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))
mod <- stan_BH
pops <- c("Lostine","Tucannon","SFSMain","SFSEast","Pahsimeroi","SREastFork")
pHOS <- extract(mod,"pHOS")[[1]]
for(i in pops)
{
  bb <- "blue4"
  plot(srchin$brood.yr[srchin$pop==i], srchin$pHOS[srchin$pop==i], 
       type="b", pch=16, cex=1.2, cex.axis=1.2, las=1, ylim=c(0,1), xlab="", ylab="")
  lines(srchin$brood.yr[stan_dat$which_H][srchin$pop[stan_dat$which_H]==i],
        apply(pHOS[,srchin$pop[stan_dat$which_H]==i], 2, quantile, 0.5), col="blue", lwd=2)
  bb <- col2rgb(bb)
  bb <- rgb(bb[1], bb[2], bb[3], maxColorValue=255, alpha=255*0.3)
  polygon(c(srchin$brood.yr[stan_dat$which_H][srchin$pop[stan_dat$which_H]==i], 
            rev(srchin$brood.yr[stan_dat$which_H][srchin$pop[stan_dat$which_H]==i])),
          c(apply(pHOS[,srchin$pop[stan_dat$which_H]==i], 2, quantile, 0.025), 
            rev(apply(pHOS[,srchin$pop[stan_dat$which_H]==i], 2, quantile, 0.975))),
          col = bb, border = NA)
  #   lines(srchin$brood.yr[srchin$pop==i],colMeans(S_tot[,srchin$pop==i]), col="blue", lwd=2)
  mtext(i, side=3, line=0.5, cex=1.2)
  if(i %in% pops[4:6]) mtext("Year", side=1, line=3, cex=1.2)
  if(i %in% pops[c(1,4)]) mtext(bquote(p[HOS]), side=2, line=3.1, cex=1.2)
}
rm(list=c("mod","pops","pHOS","bb"))
# dev.off()


# Observed spawner age proportions vs. estimated (states)
dev.new(width=10,height=10)
mod <- stan_BH
q <- matrix(get_posterior_mean(mod,"q")[,4], ncol=3, byrow=T)
q_obs <- as.matrix(srchin[,c("n3","n4","n5")])
q_obs[rowSums(q_obs)==0,] <- NA
q_obs <- sweep(q_obs, 1, rowSums(q_obs), "/")
plot(q, q_obs, pch="", xlim=c(0,1), ylim=c(0,1), cex.axis=1.2, cex.lab=1.5, cex.main=1.5, las=1, 
     xlab="Estimated", ylab="Observed", main="Spawner age distribution")
points(q[,1], q_obs[,1], pch=1, col="black", cex=1.2)
points(q[,2], q_obs[,2], pch=1, col="blue", cex=1.2)
points(q[,3], q_obs[,3], pch=1, col="green", cex=1.2)
abline(0,1)
legend("topleft", pch=1, col=c("black","blue","green"), legend = c("age 3", "age 4" ,"age 5"), pt.cex=1.2, cex=1.5)
rm(list=c("mod","q","q_obs"))


# Time series of observed and estimated spawner age proportions, panel for each of 6 pops
# png(filename="spawner_age.png", width=11, height=7, units="in", res=200, type="cairo-png")
dev.new(width=11,height=7)
par(mfrow=c(2,3), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))
mod <- stan_BH
q <- extract(mod,"q")[[1]]
q_obs <- as.matrix(srchin[,c("n3","n4","n5")])
q_obs[rowSums(q_obs)==0,] <- NA
q_obs <- sweep(q_obs, 1, rowSums(q_obs), "/")
pops <- c("Lostine","MarshCR","Tucannon","SFSEast","Pahsimeroi","SREastFork")
for(i in pops)
{
  og <- "orangered3"
  bb <- "blue4"
  gg <- "darkgray"
  plot(srchin$brood.yr[srchin$pop==i], q_obs[srchin$pop==i,1], 
       type="b", pch=16, col=og, cex=1.2, cex.axis=1.2, las=1, ylim=c(0,1), xlab="", ylab="")
  points(srchin$brood.yr[srchin$pop==i], q_obs[srchin$pop==i,2], 
         type="b", pch=16, col=bb, cex=1.2)
  points(srchin$brood.yr[srchin$pop==i], q_obs[srchin$pop==i,3], 
         type="b", pch=16, col="black", cex=1.2)
  lines(srchin$brood.yr[srchin$pop==i],apply(q[,srchin$pop==i,1], 2, quantile, 0.5), col=og, lwd=2)
  lines(srchin$brood.yr[srchin$pop==i],apply(q[,srchin$pop==i,2], 2, quantile, 0.5), col=bb, lwd=2)
  lines(srchin$brood.yr[srchin$pop==i],apply(q[,srchin$pop==i,3], 2, quantile, 0.5), col="black", lwd=2)
  og <- col2rgb(og)
  og <- rgb(og[1], og[2], og[3], maxColorValue=255, alpha=255*0.3)
  bb <- col2rgb(bb)
  bb <- rgb(bb[1], bb[2], bb[3], maxColorValue=255, alpha=255*0.3)
  gg <- col2rgb(gg)
  gg <- rgb(gg[1], gg[2], gg[3], maxColorValue=255, alpha=255*0.6)
  polygon(c(srchin$brood.yr[srchin$pop==i], rev(srchin$brood.yr[srchin$pop==i])), 
          c(apply(q[,srchin$pop==i,1], 2, quantile, 0.025), 
            rev(apply(q[,srchin$pop==i,1], 2, quantile, 0.955))),
          col = og, border = NA)
  polygon(c(srchin$brood.yr[srchin$pop==i], rev(srchin$brood.yr[srchin$pop==i])), 
          c(apply(q[,srchin$pop==i,2], 2, quantile, 0.025), 
            rev(apply(q[,srchin$pop==i,2], 2, quantile, 0.955))),
          col = bb, border = NA)
  polygon(c(srchin$brood.yr[srchin$pop==i], rev(srchin$brood.yr[srchin$pop==i])), 
          c(apply(q[,srchin$pop==i,3], 2, quantile, 0.025), 
            rev(apply(q[,srchin$pop==i,3], 2, quantile, 0.955))),
          col = gg, border = NA)
  mtext(i, side=3, line=0.5, cex=1.2)
  if(i %in% pops[4:6]) mtext("Year", side=1, line=3, cex=1.2)
  if(i %in% pops[c(1,4)]) mtext("Proportion at age", side=2, line=3.1, cex=1.2)
}
rm(list=c("mod","pops","q","q_obs","bb","og","gg"))
# dev.off()


# Estimated spawner-recruit curves (for each pop, and ESU-level)
dev.new(width=10,height=10)
par(mar=c(5.1,5.1,1,1))
mod <- stan_BH
LG <- function(a_W, a_H, Rmax_W, Rmax_H, S_W, S_H) 
{
  R_W <- a_W*S_W/(1 + a_W*S_W/Rmax_W + a_H*S_H/Rmax_H);
  R_H <- a_H*S_H/(1 + a_W*S_W/Rmax_W + a_H*S_H/Rmax_H);
  R <- R_W + R_H;
  return(R);
}
S_W <- matrix(seq(0,quantile(stan_dat$S_tot_obs/stan_dat$A, 0.9),length=100),
              nrow=sum(stan_BH@sim$n_save - stan_BH@sim$warmup2), ncol=100, byrow=T)
R_ESU <- LG(a_W = as.vector(exp(extract(mod,"mu_log_a_W")[[1]])), a_H = 0, 
            Rmax_W = as.vector(exp(extract(mod,"mu_log_Rmax_W")[[1]])), Rmax_H = 1, 
            S_W = S_W, S_H = 0)
plot(S_W[1,], apply(R_ESU,2,median), type="l", col="blue", lwd=3, las=1, 
     xaxs="i", yaxs="i", ylim=c(0,30), cex.lab=1.5, cex.axis=1.2,
     xlab=bquote("Spawner density (ha"^-1*")"), ylab=bquote("Recruit density (ha"^-1*")"))
for(j in 1:length(levels(srchin$pop)))
{
  R <- LG(a_W = as.vector(extract(mod,"a_W")[[1]][,j]), a_H = 0, 
          Rmax_W = as.vector(extract(mod,"Rmax_W")[[1]][,j]), Rmax_H = 1, 
          S_W = S_W, S_H = 0)
  lines(S_W[1,], apply(R,2,median), col="darkgray", lwd=1.5)
}
bb <- col2rgb("darkblue")
bb <- rgb(bb[1], bb[2], bb[3], maxColorValue = 255, alpha = 255*0.3)
polygon(c(S_W[1,], rev(S_W[1,])), 
        c(apply(R_ESU,2,quantile,0.025), rev(apply(R_ESU,2,quantile,0.975))),
        col=bb, border=NA)
rm(list=c("S_W","R_ESU","R","LG","mod","bb"))


# Spawner-recruit curve variance decomposition:
# estimation error, spatial (among pops), temporal (among years), process, and obs error
# png(filename="state-space_SR.png", width=7, height=7, units="in", res=200, type="cairo-png")
dev.new(width=10,height=10)
par(mar=c(5.1,5.1,1,1))
LG <- function(a_W, a_H, Rmax_W, Rmax_H, S_W, S_H) 
{
  R_W <- a_W*S_W/(1 + a_W*S_W/Rmax_W + a_H*S_H/Rmax_H);
  R_H <- a_H*S_H/(1 + a_W*S_W/Rmax_W + a_H*S_H/Rmax_H);
  R <- R_W + R_H;
  return(R);
}
mod <- stan_BH
mu_log_a_W <- as.vector(extract(mod,"mu_log_a_W")[[1]])
sigma_log_a <- as.vector(extract(mod,"sigma_log_a")[[1]])
mu_log_Rmax_W <- as.vector(extract(mod,"mu_log_Rmax_W")[[1]])
sigma_log_Rmax <- as.vector(extract(mod,"sigma_log_Rmax")[[1]])
sigma_log_phi <- as.vector(extract(mod,"sigma_log_phi")[[1]])
sigma_proc <- as.vector(extract(mod,"sigma_proc")[[1]])
sigma_obs <- as.vector(extract(mod,"sigma_obs")[[1]])
S_W <- matrix(seq(0,quantile(stan_dat$S_tot_obs/stan_dat$A, 0.9),length=100),
              nrow=sum(mod@sim$n_save - mod@sim$warmup2)*50, ncol=100, byrow=T)
R_ESU <- LG(a_W = exp(mu_log_a_W), a_H = 0, 
            Rmax_W = exp(mu_log_Rmax_W), Rmax_H = 1, 
            S_W = S_W, S_H = 0)
R_pop <- LG(a_W = rlnorm(nrow(S_W), mu_log_a_W, sigma_log_a), a_H = 0, 
            Rmax_W = rlnorm(nrow(S_W), mu_log_Rmax_W, sigma_log_Rmax), Rmax_H = 1, 
            S_W = S_W, S_H = 0)
R_year <- R_pop*rlnorm(nrow(S_W), 0, sigma_log_phi)
R_proc <- R_year*rlnorm(nrow(S_W), 0, sigma_proc)
R_obs <- R_proc*rlnorm(nrow(S_W), 0, sigma_obs)
bb <- "blue4"
plot(S_W[1,], apply(R_obs,2,median), type="l", lwd=3, col=bb, las=1,
     cex.lab=1.5, cex.axis=1.2, xaxs="i", yaxs="i", 
     ylim=c(0,max(apply(R_obs,2,quantile,0.975))),
     xlab=bquote("Spawner density (ha"^-1*")"), ylab=bquote("Recruit density (ha"^-1*")"))
bb <- col2rgb(bb)
bb <- rgb(bb[1], bb[2], bb[3], maxColorValue = 255, alpha = 255*0.15)
polygon(c(S_W[1,], rev(S_W[1,])), 
        c(apply(R_ESU,2,quantile,0.025), rev(apply(R_ESU,2,quantile,0.975))), col=bb, border=NA)
polygon(c(S_W[1,], rev(S_W[1,])), 
        c(apply(R_pop,2,quantile,0.025), rev(apply(R_pop,2,quantile,0.975))), col=bb, border=NA)
polygon(c(S_W[1,], rev(S_W[1,])), 
        c(apply(R_year,2,quantile,0.025), rev(apply(R_year,2,quantile,0.975))), col=bb, border=NA)
polygon(c(S_W[1,], rev(S_W[1,])), 
        c(apply(R_proc,2,quantile,0.025), rev(apply(R_proc,2,quantile,0.975))), col=bb, border=NA)
polygon(c(S_W[1,], rev(S_W[1,])), 
        c(apply(R_obs,2,quantile,0.025), rev(apply(R_obs,2,quantile,0.975))), col=bb, border=NA)
rm(list=c("mod","mu_log_a_W","sigma_log_a","mu_log_Rmax_W","sigma_log_Rmax","sigma_log_phi",
          "sigma_proc","sigma_obs","S_W","R_ESU","R_pop","R_year","R_proc","R_obs","bb"))
# dev.off()


# TRADITIONAL REGRESSION-BASED spawner-recruit curve variance decomposition:
# estimation error, spatial (among pops), temporal (among years), and residual error
# png(filename="regression_SR.png", width=7, height=7, units="in", res=200, type="cairo-png")
dev.new(width=10,height=10)
par(mar=c(5.1,5.1,1,1))
LG <- function(a_W, a_H, Rmax_W, Rmax_H, S_W, S_H) 
{
  R_W <- a_W*S_W/(1 + a_W*S_W/Rmax_W + a_H*S_H/Rmax_H);
  R_H <- a_H*S_H/(1 + a_W*S_W/Rmax_W + a_H*S_H/Rmax_H);
  R <- R_W + R_H;
  return(R);
}
mod <- stan_BH_reg
mu_log_a_W <- as.vector(extract(mod,"mu_log_a_W")[[1]])
sigma_log_a <- as.vector(extract(mod,"sigma_log_a")[[1]])
mu_log_Rmax_W <- as.vector(extract(mod,"mu_log_Rmax_W")[[1]])
sigma_log_Rmax <- as.vector(extract(mod,"sigma_log_Rmax")[[1]])
sigma_log_phi <- as.vector(extract(mod,"sigma_log_phi")[[1]])
sigma <- as.vector(extract(mod,"sigma")[[1]])
S_W <- matrix(seq(0,quantile((stan_reg_dat$S_W + stan_reg_dat$S_H)/stan_reg_dat$A, 0.9),length=100),
              nrow=sum(mod@sim$n_save - mod@sim$warmup2)*50, ncol=100, byrow=T)
R_ESU <- LG(a_W = exp(mu_log_a_W), a_H = 0, 
            Rmax_W = exp(mu_log_Rmax_W), Rmax_H = 1, 
            S_W = S_W, S_H = 0)
R_pop <- LG(a_W = rlnorm(nrow(S_W), mu_log_a_W, sigma_log_a), a_H = 0, 
            Rmax_W = rlnorm(nrow(S_W), mu_log_Rmax_W, sigma_log_Rmax), Rmax_H = 1, 
            S_W = S_W, S_H = 0)
R_year <- R_pop*rlnorm(nrow(S_W), 0, sigma_log_phi)
R_resid <- R_year*rlnorm(nrow(S_W), 0, sigma)
bb <- "orangered3"
plot(S_W[1,], apply(R_resid,2,median), type="l", lwd=3, col=bb, las=1,
     cex.lab=1.5, cex.axis=1.2, xaxs="i", yaxs="i", 
     ylim=c(0,216),
     xlab=bquote("Spawner density (ha"^-1*")"), ylab=bquote("Recruit density (ha"^-1*")"))
bb <- col2rgb(bb)
bb <- rgb(bb[1], bb[2], bb[3], maxColorValue = 255, alpha = 255*0.15)
polygon(c(S_W[1,], rev(S_W[1,])), 
        c(apply(R_ESU,2,quantile,0.025), rev(apply(R_ESU,2,quantile,0.975))), col=bb, border=NA)
polygon(c(S_W[1,], rev(S_W[1,])), 
        c(apply(R_pop,2,quantile,0.025), rev(apply(R_pop,2,quantile,0.975))), col=bb, border=NA)
polygon(c(S_W[1,], rev(S_W[1,])), 
        c(apply(R_year,2,quantile,0.025), rev(apply(R_year,2,quantile,0.975))), col=bb, border=NA)
polygon(c(S_W[1,], rev(S_W[1,])), 
        c(apply(R_resid,2,quantile,0.025), rev(apply(R_resid,2,quantile,0.975))), col=bb, border=NA)
rm(list=c("mod","mu_log_a_W","sigma_log_a","mu_log_Rmax_W","sigma_log_Rmax","sigma_log_phi",
          "sigma","S_W","R_ESU","R_pop","R_year","R_resid","bb"))
# dev.off()


# Comparison of "spawner-recruit" regression vs. integrated model S-R curves
# png(filename="regression_vs_SS_SR.png", width=7, height=7, units="in", res=200, type="cairo-png")
dev.new(width=10,height=10)
par(mar=c(5.1,5.1,1,1))
LG <- function(a_W, a_H, Rmax_W, Rmax_H, S_W, S_H) 
{
  R_W <- a_W*S_W/(1 + a_W*S_W/Rmax_W + a_H*S_H/Rmax_H);
  R_H <- a_H*S_H/(1 + a_W*S_W/Rmax_W + a_H*S_H/Rmax_H);
  R <- R_W + R_H;
  return(R);
}
mod <- stan_BH
mu_log_a_W <- as.vector(extract(mod,"mu_log_a_W")[[1]])
mu_log_Rmax_W <- as.vector(extract(mod,"mu_log_Rmax_W")[[1]])
S_W <- matrix(seq(0,quantile((stan_dat$S_tot_obs)/stan_dat$A, 0.9),length=100),
              nrow=sum(mod@sim$n_save - mod@sim$warmup2)*50, ncol=100, byrow=T)
R_ESU <- LG(a_W = exp(mu_log_a_W), a_H = 0, 
            Rmax_W = exp(mu_log_Rmax_W), Rmax_H = 1, 
            S_W = S_W, S_H = 0)
bb <- "blue4"
plot(S_W[1,], apply(R_ESU,2,median), type="l", lwd=3, col=bb, las=1,
     cex.lab=1.5, cex.axis=1.2, xaxs="i", yaxs="i", 
     ylim=c(0,max(apply(R_ESU,2,quantile,0.975))),
     xlab=bquote("Spawner density (ha"^-1*")"), ylab=bquote("Recruit density (ha"^-1*")"))
bb <- col2rgb(bb)
bb <- rgb(bb[1], bb[2], bb[3], maxColorValue = 255, alpha = 255*0.15)
polygon(c(S_W[1,], rev(S_W[1,])), 
        c(apply(R_ESU,2,quantile,0.025), rev(apply(R_ESU,2,quantile,0.975))), col=bb, border=NA)
rm(list=c("mod","mu_log_a_W","mu_log_Rmax_W","R_ESU","bb"))
mod <- stan_BH_reg
mu_log_a_W <- as.vector(extract(mod,"mu_log_a_W")[[1]])
mu_log_Rmax_W <- as.vector(extract(mod,"mu_log_Rmax_W")[[1]])
R_ESU <- LG(a_W = exp(mu_log_a_W), a_H = 0, 
            Rmax_W = exp(mu_log_Rmax_W), Rmax_H = 1, 
            S_W = S_W, S_H = 0)
bb <- "orangered3"
lines(S_W[1,], apply(R_ESU,2,median), type="l", lwd=3, col=bb)
bb <- col2rgb(bb)
bb <- rgb(bb[1], bb[2], bb[3], maxColorValue = 255, alpha = 255*0.15)
polygon(c(S_W[1,], rev(S_W[1,])), 
        c(apply(R_ESU,2,quantile,0.025), rev(apply(R_ESU,2,quantile,0.975))), col=bb, border=NA)
rm(list=c("mod","mu_log_a_W","mu_log_Rmax_W","S_W","R_ESU","bb"))
# dev.off()


#-------------------------------------------------------
# Diagnostic plots
# (Require monitoring log_phi_z, log_R_tot_z, alr_p_z) 
#-------------------------------------------------------

# Residuals of log(phi): normality and autocorrelation
dev.new(width=12,height=7)
par(mfrow=c(1,2))
qqnorm(stan_mean(stan_BH,"log_phi_z"))
qqline(stan_mean(stan_BH,"log_phi_z"))
acf(stan_mean(stan_BH,"log_phi_z"))

# Recruitment process error residuals, log(R_tot) - log(R_tot_hat): 
# normality and autocorrelation w/in pops,
# cross-correlation b/w pops, boxplots by pop
dev.new()
qqnorm(stan_mean(stan_BH,"log_R_tot_z"))
qqline(stan_mean(stan_BH,"log_R_tot_z"))

dev.new(width=12,height=7)
par(mfrow=c(4,6), mar = c(4.1,4.1,4.1,1.1), oma=c(0,0,3,0))
for(i in levels(srchin$pop))
  acf(stan_mean(stan_BH,"log_R_tot_z")[srchin$pop==i], main=i)
mtext("log(R_tot) - log(R_tot_hat)", outer = TRUE)

dev.new()
bwplot(srchin$pop ~ stan_mean(stan_BH,"log_R_tot_z"))

# Age structure process error residuals, alr(p[i,age]) - mu_alr_p[age]:
# normality and autocorrelation w/in pops,
# cross-correlation b/w pops,
# cross-correlation b/w ages,
# boxplots by pop
alr_p_z <- matrix(stan_mean(stan_BH,"alr_p_z"), ncol = 2, byrow = TRUE)

dev.new(width=12,height=7)
par(mfrow=c(1,2))
qqnorm(alr_p_z[,1])
qqline(alr_p_z[,1])
qqnorm(alr_p_z[,2])
qqline(alr_p_z[,2])

dev.new(width=12,height=7)
par(mfrow=c(4,6), mar = c(4.1,4.1,4.1,1.1), oma=c(0,0,3,0))
for(i in levels(srchin$pop))
  acf(alr_p_z[srchin$pop==i,1], main=i)
mtext("log(p3/p5)", outer = TRUE)

dev.new(width=12,height=7)
par(mfrow=c(4,6), mar = c(4.1,4.1,4.1,1.1), oma=c(0,0,3,0))
for(i in levels(srchin$pop))
  acf(alr_p_z[srchin$pop==i,2], main=i)
mtext("log(p4/p5)", outer = TRUE)

dev.new(width=14,height=10)
par(mfrow=c(4,6), mar = c(4.1,4.1,4.1,1.1))
for(i in levels(srchin$pop))
  plot(alr_p_z[srchin$pop==i,1], alr_p_z[srchin$pop==i,2], xlab="log(p3/p5)", ylab="log(p4/p5)", main=i)

dev.new()
bwplot(srchin$pop ~ alr_p_z[,1])

dev.new()
bwplot(srchin$pop ~ alr_p_z[,2])

rm(alr_p_z)

# Observation error, log(S_tot_obs) - log(S_tot):
# normality and autocorrelation w/in pops,
# boxplots by pop
obs_err <- log(stan_dat$S_tot_obs) - colMeans(log(extract1(stan_BH,"S_tot")))

dev.new()
qqnorm(obs_err)
qqline(obs_err)

dev.new(width=12,height=7)
par(mfrow=c(4,6), mar = c(4.1,4.1,4.1,1.1), oma=c(0,0,3,0))
for(i in levels(srchin$pop))
  acf(obs_err[srchin$pop==i], main=i)
mtext("obs error", outer = TRUE)

dev.new()
bwplot(srchin$pop ~ obs_err)

rm(obs_err)





