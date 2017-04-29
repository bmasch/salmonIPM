options(device=windows)
library(salmonIPM)

#===========================================================================
# DATA
#===========================================================================

# Load data
fish_data <- read.table(file.path("~", "SalmonIPM", "IPM_PVA", "fish_data.txt"), sep = "\t", header = T)
fish_data <- fish_data[order(fish_data$code, fish_data$year),]

# Impute one NA value of S_tot_obs in Chamberlain 1986
fish_data$S_tot_obs[fish_data$pop == "Chamberlain" & fish_data$year == 1986] <- 
  mean(fish_data$S_tot_obs[fish_data$pop == "Chamberlain"][1:5], na.rm = T)

# Change area to 1 for all pops (units of Rmax will be spawners, not spawners/ha)
# fish_data$A <- 1

# Pad data with years through max_year
max_year <- max(fish_data$year) + 50
year_aug <- sapply(tapply(fish_data$year, fish_data$pop, max), function(x) (x + 1):max_year)
pop_aug <- rep(names(year_aug), sapply(year_aug, length))
code_aug <- fish_data$code[match(pop_aug, fish_data$pop)]
MPG_aug <- fish_data$MPG[match(pop_aug, fish_data$pop)]
A_aug <- rep(tapply(fish_data$A, fish_data$pop, mean), times = sapply(year_aug, length))
fish_data_aug <- data.frame(pop = pop_aug, code = code_aug, MPG = MPG_aug, A = A_aug,
                            year = unlist(year_aug), type = "future", fit_p_HOS = 0, 
                            S_tot_obs = NA, n_age3_obs = 0, n_age4_obs = 0, n_age5_obs = 0,
                            n_W_obs = 0, n_H_obs = 0, p_HOS = 0, B_take_obs = 0, F_rate = 0,
                            row.names = NULL)
fish_data_aug <- rbind(cbind(type = "past", fish_data[,setdiff(names(fish_data_aug), "type")])[,names(fish_data_aug)], 
                       fish_data_aug)
fish_data_aug <- fish_data_aug[order(fish_data_aug$code, fish_data_aug$year),]
row.names(fish_data_aug) <- NULL


#===========================================================================
# FIT MODELS
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

print(PVA_IPM_pp, pars = c("phi","p_HOS","B_rate_all","q"), include = FALSE)
launch_shinystan(PVA_IPM_pp)









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
# FIGURES
#===========================================================================

# Comparison of S-R curves and parameters under RR and IPM models
dev.new(width = 15, height = 5)
par(mfrow = c(1,3))
BH <- function(a, b, S) 
{
  a*S/(1 + b*S)
}

S <- matrix(seq(1, quantile(fish_data$S_tot_obs, 0.9, na.rm = T), length = 100),
            nrow = sum(PVA_RR_pp@sim$n_save - PVA_RR_pp@sim$warmup2)*50, ncol = 100, byrow = T)

# S-R curves
mu_log_a <- as.vector(extract1(PVA_IPM_pp,"mu_log_a"))
mu_log_b <- as.vector(extract1(PVA_IPM_pp,"mu_log_b"))
R_ESU <- BH(a = exp(mu_log_a), b = exp(mu_log_b), S = S)
bb <- "blue4"

plot(S[1,], apply(R_ESU, 2, median), type = "l", lwd=3, col = bb, las = 1,
     cex.lab = 1.5, cex.axis = 1.2, xaxs = "i", yaxs = "i",
     ylim = c(0, max(apply(R_ESU, 2, quantile, 0.975))),
     xlab = "Spawners", ylab = "Recruits")
bb <- col2rgb(bb)
bb <- rgb(bb[1], bb[2], bb[3], maxColorValue = 255, alpha = 255*0.15)
polygon(c(S[1,], rev(S[1,])), c(apply(R_ESU, 2, quantile, 0.025), rev(apply(R_ESU, 2, quantile, 0.975))), col = bb, border = NA)

mu_log_a <- as.vector(extract1(PVA_RR_pp,"mu_log_a"))
mu_log_b <- as.vector(extract1(PVA_RR_pp,"mu_log_b"))
R_ESU <- BH(a = exp(mu_log_a), b = exp(mu_log_b), S = S)
bb <- "orangered3"
lines(S[1,], apply(R_ESU, 2, median), type = "l", lwd=3, col = bb)
bb <- col2rgb(bb)
bb <- rgb(bb[1], bb[2], bb[3], maxColorValue = 255, alpha = 255*0.15)
polygon(c(S[1,], rev(S[1,])), c(apply(R_ESU, 2, quantile, 0.025), rev(apply(R_ESU, 2, quantile, 0.975))), col = bb, border = NA)

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
plot(dd_IPM_ESU$x, dd_IPM_ESU$y, type = "l", lwd = 3, col = bb, las = 1, cex.lab = 1.5, cex.axis = 1.2,
     xlab = "log(a)", ylab = "Probability density", xaxs = "i",
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

# Posterior densities of log(Rmax)
dd_IPM_ESU <- density(extract1(PVA_IPM_pp,"mu_log_a") - extract1(PVA_IPM_pp,"mu_log_b"))
dd_IPM_pop <- vector("list", length(levels(fish_data$pop)))
for(i in 1:length(dd_IPM_pop))
  dd_IPM_pop[[i]] <- density(log(extract1(PVA_IPM_pp,"a")[,i]) - log(extract1(PVA_IPM_pp,"b")[,i]))
dd_RR_ESU <- density(extract1(PVA_RR_pp,"mu_log_a") - extract1(PVA_RR_pp,"mu_log_b"))
dd_RR_pop <- vector("list", length(levels(fish_data$pop)))
for(i in 1:length(dd_RR_pop))
  dd_RR_pop[[i]] <- density(log(extract1(PVA_RR_pp,"a")[,i]) - log(extract1(PVA_RR_pp,"b")[,i]))

bb <- "blue4"
plot(dd_IPM_ESU$x, dd_IPM_ESU$y, type = "l", lwd = 3, col = bb, las = 1, cex.lab = 1.5, cex.axis = 1.2,
     xlab = "log(Rmax)", ylab = "Probability density", xaxs = "i",
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

rm(list=c("mu_log_a","sigma_log_a","mu_log_b","sigma_log_b","S","R_ESU",
          "bb","dd_IPM_ESU","dd_RR_ESU","dd_IPM_pop","dd_RR_pop"))



# Spawner-recruit curve variance decomposition under RR and IPM models

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

rm(list=c("mu_log_a","sigma_log_a","mu_log_b","sigma_log_b","sigma_log_phi","sigma",
          "mu_sigma_proc","sigma_obs","S","R_ESU","R_pop","R_year","R_proc","R_obs","R_resid","bb"))




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




