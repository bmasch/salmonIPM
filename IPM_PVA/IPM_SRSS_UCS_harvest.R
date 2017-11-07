setwd(file.path("~","salmonIPM","IPM_PVA"))
options(device=windows)
library(salmonIPM)
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

# Create "forward simulation" data
F_rate_fwd <- seq(0, 0.3, by = 0.05)
N_future_years <- 50
max_year <- max(fish_data$year) + N_future_years
year_fwd <- sapply(tapply(fish_data$year, fish_data$pop, max), function(x) (x + 1):max_year)
pop_fwd <- rep(names(year_fwd), sapply(year_fwd, length))
A_fwd <- rep(tapply(fish_data$A, fish_data$pop, mean), times = sapply(year_fwd, length))
year_fwd = unlist(year_fwd)
fish_data_fwd <- data.frame(pop = rep(pop_fwd, length(F_rate_fwd)), 
                            year = rep(year_fwd, length(F_rate_fwd)), 
                            A = rep(A_fwd, length(F_rate_fwd)), 
                            F_rate = rep(F_rate_fwd, each = length(year_fwd)), 
                            B_rate = 0, p_HOS = 0)
row.names(fish_data_fwd) <- NULL

# Create dummy covariate for intervention analysis: 
# pre/post-1970 (centered, not scaled, based on years 1:N_year)
# years (N_year + 1):N_year_all are given the post-1970 value
years <- sort(unique(c(fish_data$year, fish_data_fwd$year)))
step_1970 <- data.frame(step_1970 = as.numeric(years >= 1970),
                            row.names = years)
step_1970$step_1970 <- step_1970$step_1970 - 
  mean(step_1970[as.character(min(fish_data$year):max(fish_data$year)),"step_1970"])


#===========================================================================
# HARVEST SIMULATIONS
# Vary future exploitation rate (fixed-rate policy) and look at
# impact on quasi-extinction risk under multi-pop IPM
#===========================================================================

#--------------------------------------------------------------------------
# Model without any covariates
#--------------------------------------------------------------------------

PVA_F <- salmonIPM(fish_data = fish_data, fish_data_fwd = fish_data_fwd, model = "IPM", pool_pops = TRUE, 
                   pars = c("mu_log_a","sigma_log_a","a",
                            "mu_log_Rmax","sigma_log_Rmax","Rmax","rho_log_aRmax",
                            "beta_log_phi","sigma_log_phi","rho_log_phi","phi",
                            "mu_p","sigma_gamma","R_gamma","gamma",
                            "sigma_alr_p","R_alr_p","p","p_fwd",
                            "p_HOS","sigma_proc","sigma_obs",
                            "S_tot","S_tot_fwd","R_tot","R_tot_fwd","q","q_fwd"),
                   chains = 3, iter = 1500, warmup = 500,
                   control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

print(PVA_F, pars = c("a","Rmax","phi","p_HOS","B_rate_all","q","q_fwd","gamma",
                      "p","p_fwd","S_tot","S_tot_fwd","R_tot","R_tot_fwd"), 
      include = FALSE, use_cache = FALSE)

launch_shinystan(as.shinystan(PVA_F,
                              pars = c("mu_log_a","sigma_log_a","a",
                                       "mu_log_Rmax","sigma_log_Rmax","Rmax","rho_log_aRmax",
                                       "beta_log_phi","sigma_log_phi","rho_log_phi","phi",
                                       "mu_p","sigma_gamma","R_gamma","gamma",
                                       "sigma_alr_p","R_alr_p","sigma_proc","sigma_obs")))

# Convert to array and overwrite stanfit object to save space
### (ARE YOU SURE YOU'RE READY TO DO THIS??)
PVA_F <- as.array(PVA_F)
gc()


#--------------------------------------------------------------------------
# Model with a step change in 1970
#--------------------------------------------------------------------------

PVA_F_1970 <- salmonIPM(fish_data = fish_data, fish_data_fwd = fish_data_fwd, env_data = step_1970,
                        model = "IPM", pool_pops = TRUE, 
                        pars = c("mu_log_a","sigma_log_a","a",
                                 "mu_log_Rmax","sigma_log_Rmax","Rmax","rho_log_aRmax",
                                 "beta_log_phi","sigma_log_phi","rho_log_phi","phi",
                                 "mu_p","sigma_gamma","R_gamma","gamma",
                                 "sigma_alr_p","R_alr_p","p","p_fwd",
                                 "p_HOS","sigma_proc","sigma_obs",
                                 "S_tot","S_tot_fwd","R_tot","R_tot_fwd","q","q_fwd"),
                        chains = 3, iter = 1500, warmup = 500,
                        control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

print(PVA_F_1970, pars = c("a","Rmax","phi","p_HOS","B_rate_all","q","q_fwd","gamma",
                      "p","p_fwd","S_tot","S_tot_fwd","R_tot","R_tot_fwd"), 
      include = FALSE, use_cache = FALSE)

launch_shinystan(as.shinystan(PVA_F_1970, 
                              pars = c("mu_log_a","sigma_log_a","a",
                                       "mu_log_Rmax","sigma_log_Rmax","Rmax","rho_log_aRmax",
                                       "beta_log_phi","sigma_log_phi","rho_log_phi","phi",
                                       "mu_p","sigma_gamma","R_gamma","gamma",
                                       "sigma_alr_p","R_alr_p","sigma_proc","sigma_obs")))

# Convert to array and delete stanfit object to save space
### (ARE YOU SURE YOU'RE READY TO DO THIS??)
PVA_F_1970 <- as.array(PVA_F_1970)
gc()
                 

#===========================================================================
# FIGURES
#===========================================================================

#------------------------------------------------------------------------------------
# Sample paths of brood-year log productivity anomalies (log phi)
#------------------------------------------------------------------------------------

dev.new(width = 14, height = 7)
png(filename="log_phi_sample_paths.png", width=14, height=7, units="in", res=200, type="cairo-png")
layout(matrix(c(1,2), nrow = 1), widths = c(14,2))

mod <- "PVA_F_1970"
type <- ifelse(years %in% fish_data$year, "past", "future")
log_phi <- log(eval(as.name(mod))[,,substring(dimnames(eval(as.name(mod)))[[3]],1,4)=="phi["])
log_phi <- matrix(log_phi, nrow = prod(dim(log_phi)[1:2]))
beta_log_phi <- as.vector(eval(as.name(mod))[,,grep("beta_log_phi", dimnames(eval(as.name(mod)))[[3]])])
quantiles <- quantile(rowMeans(log_phi[,type == "future"]), c(1/3, 2/3))
gmb <- ifelse(rowMeans(log_phi[,type == "future"]) >= quantiles[1],
              ifelse(rowMeans(log_phi[,type == "future"]) >= quantiles[2], "good", "medium"),
              "bad")

c1 <- "blue4"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.5)
c2 <- "darkgray"
c2t <- col2rgb(c2)
c2t <- rgb(c2t[1], c2t[2], c2t[3], maxColorValue = 255, alpha = 255*0.5)
alpha <- c(0.7, 0.8, 0.5)
cgmb <- lapply(c("green", "yellow2", "red"), col2rgb)
cgmb <- sapply(c(good = 1, medium = 2, bad = 3), 
               function(i) rgb(cgmb[[i]][1], cgmb[[i]][2], cgmb[[i]][3], maxColorValue = 255, alpha = 255*alpha[i]))

indx1 <- c(which.max(rowMeans(log_phi[,type=="future"])), which.min(rowMeans(log_phi[,type=="future"])))
indx2 <- sample(nrow(log_phi),100)

# layer 1
par(mar = c(5.1, 4.5, 4.1, 0.2))
plot(years, log_phi[indx1[1],], type = "l", col = c1t, las = 1, cex.lab = 1.8, cex.axis = 1.5, cex.main = 1.8,
     xaxt = "n", xaxs = "i", yaxs = "i", ylim = c(-5,5), #range(log_phi[union(indx1,indx2),]), 
     xlab = "Brood year", ylab = "Productivity anomaly", 
     main = ifelse(mod == "PVA_F", "Constant baseline", "Step change 1970"))
abline(v = max(years[type == "past"]), lty = 2)
if(mod == "PVA_F") abline(h = 0, col = c2, lwd = 3)
if(mod  ==  "PVA_F_1970")
{
  mu_log_phi <- as.matrix(beta_log_phi) %*% t(step_1970)
  lines(years[years < 1970], colMeans(mu_log_phi[,years < 1970]), col = c2, lwd = 3)
  polygon(c(years[years < 1970], rev(years[years < 1970])),
          c(apply(mu_log_phi[,years < 1970], 2, quantile, 0.025),
            rev(apply(mu_log_phi[,years < 1970], 2, quantile, 0.975))), col = c2t, border = NA)
  lines(years[years >= 1970], colMeans(mu_log_phi[,years >= 1970]), col = c2, lwd = 3)
  polygon(c(years[years >= 1970], rev(years[years >= 1970])),
          c(apply(mu_log_phi[,years >= 1970], 2, quantile, 0.025),
            rev(apply(mu_log_phi[,years >= 1970], 2, quantile, 0.975))), col = c2t, border = NA)
  rm(mu_log_phi)
}
axis(side = 1, at = years[years %% 10 == 0], cex.axis = 1.5)
rug(years[years %% 10 != 0], ticksize = -0.01)
rug(years[years %% 10 != 0 & years %% 5 == 0], ticksize = -0.02)

# layer 2
lines(years, log_phi[indx1[2],], col = c1t)

# layer 3
lines(years[type=="future"], log_phi[indx1[1],type=="future"], col = "green")
lines(years[type=="future"], log_phi[indx1[2],type=="future"], col = "red")
segments(x0 = max(years[type == "past"]), x1 = max(years), 
         y0 = rowMeans(log_phi[indx1,type=="future"]), col = c("green","red"))

# layer 4
for(i in indx2)
{
  lines(years[type=="past"], log_phi[i,type=="past"], col = c1t)
  lines(c(max(years[type=="past"]), years[type=="future"]),
        c(log_phi[i,sum(type=="past")], log_phi[i,type=="future"]), col = cgmb[[gmb[i]]])
}
dd <- density(rowMeans(log_phi[,type=="future"]))
par(mar = c(5.1, 0, 4.1, 0.5))
plot(dd$y, dd$x, type = "l", ylim = par("usr")[3:4], col = c1, lwd = 2,
     xlab = "", ylab = "", xlim = c(-0.01,max(dd$y)*1.05), xaxs = "i", xaxt = "n", yaxt = "n")
x1 <- which.min(abs(dd$x - quantiles[1]))
x2 <- which.min(abs(dd$x - quantiles[2]))
x3 <- length(dd$x)
polygon(c(rep(0,x1), dd$y[1:x1]), c(dd$x[x1:1], dd$x[1:x1]), col = "red", border = NA)
polygon(c(rep(0, x2 - x1), dd$y[(x1+1):x2]), c(dd$x[x2:(x1+1)], dd$x[(x1+1):x2]), col = "yellow", border = NA)
polygon(c(rep(0, x3 - x2), dd$y[(x2+1):x3]), c(dd$x[x3:(x2+1)], dd$x[(x2+1):x3]), col = "green", border = NA)
segments(x0 = 0, x1 = dd$y[which.min(abs(dd$x - quantiles[1]))], y0 = quantiles[1], col = c1)
segments(x0 = 0, x1 = dd$y[which.min(abs(dd$x - quantiles[2]))], y0 = quantiles[2], col = c1)

rm(list = c("mod","type","log_phi","beta_log_phi","c1","c1t","c2","c2t","gmb","cgmb","alpha",
            "indx1","indx2","dd","quantiles","x1","x2","x3"))
dev.off()


#------------------------------------------------------------------------------------
# Sample paths of S_tot for one population
#------------------------------------------------------------------------------------

dev.new(width = 14, height = 7)
png(filename="S_tot_sample_paths.png", width=14, height=7, units="in", res=200, type="cairo-png")
layout(matrix(c(1,2), nrow = 1), widths = c(14,2))

mod <- "PVA_F_1970"
pop <- "Lemhi"
S_tot <- eval(as.name(mod))[,,substring(dimnames(eval(as.name(mod)))[[3]],1,6)=="S_tot["]
S_tot <- S_tot[,,fish_data$pop==pop]
S_tot <- matrix(S_tot, nrow = prod(dim(S_tot)[1:2]))
S_tot_fwd <- eval(as.name(mod))[,,substring(dimnames(eval(as.name(mod)))[[3]],1,9)=="S_tot_fwd"]
S_tot_fwd <- S_tot_fwd[,,fish_data_fwd$pop==pop & fish_data_fwd$F_rate==0]
S_tot_fwd <- matrix(S_tot_fwd, nrow = prod(dim(S_tot_fwd)[1:2]))
S_tot_fwd <- cbind(S_tot[,ncol(S_tot)], S_tot_fwd)
type <- ifelse(years %in% fish_data$year, "past", "future")
log_phi <- log(eval(as.name(mod))[,,substring(dimnames(eval(as.name(mod)))[[3]],1,4)=="phi["])
log_phi <- matrix(log_phi, nrow = prod(dim(log_phi)[1:2]))
quantiles <- quantile(rowMeans(log_phi[,type == "future"]), c(1/3, 2/3))
indx <- sample(nrow(S_tot),100)
yy <- sort(unique(c(fish_data$year[fish_data$pop==pop], fish_data_fwd$year[fish_data_fwd$pop==pop])))
gmb <- ifelse(rowMeans(log_phi[,type == "future"]) >= quantiles[1],
              ifelse(rowMeans(log_phi[,type == "future"]) >= quantiles[2], "good", "medium"),
              "bad")

c1 <- "blue4"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.3)
alpha <- c(0.7, 0.8, 0.5)
cgmb <- lapply(c("green", "yellow2", "red"), col2rgb)
cgmb <- lapply(c(good = 1, medium = 2, bad = 3), 
               function(i) rgb(cgmb[[i]][1], cgmb[[i]][2], cgmb[[i]][3], maxColorValue = 255, alpha = 255*alpha[i]))

# panel 1: timeseries
par(mar = c(5.1, 4.5, 4.1, 0.2))
plot(fish_data$year[fish_data$pop==pop], fish_data$S_tot_obs[fish_data$pop==pop], pch = "",
     xlim = c(min(yy) - 0.5, max(yy)), ylim = c(1,1e5), #c(1, max(S_tot_fwd[indx,])), 
     cex.axis = 1.2, cex.lab = 1.8, cex.main = 1.8, las = 1, 
     xaxs = "i", xaxt = "n", yaxt = "n", log = "y",
     xlab = "Year", ylab = "Spawners", 
     main = paste0(pop, " (", ifelse(mod == "PVA_F", "Constant baseline", "Step change 1970"), ")"))
abline(v = max(fish_data$year[fish_data$pop==pop]), lty = 2)
axis(side = 1, at = years[years %% 10 == 0], cex.axis = 1.5)
rug(yy[yy %% 10 != 0], ticksize = -0.01)
rug(yy[yy %% 10 != 0 & yy %% 5 == 0], ticksize = -0.02)
at <- maglab(10^par("usr")[3:4], log = T)
axis(2, at$labat, cex.axis=1.5, las=1,
     labels = sapply(log10(at$labat), function(i) as.expression(bquote(10^ .(i)))))
for(i in indx)
{
  lines(fish_data$year[fish_data$pop==pop], S_tot[i,], col = c1t)
  lines(c(max(fish_data$year[fish_data$pop==pop]),
          fish_data_fwd$year[fish_data_fwd$pop==pop & fish_data_fwd$F_rate==0]), 
        S_tot_fwd[i,], col = cgmb[[gmb[i]]])
}
points(fish_data$year[fish_data$pop==pop], fish_data$S_tot_obs[fish_data$pop==pop], pch=16, cex = 1.2)

# panel 2: densities
ddg <- density(log10(S_tot_fwd[gmb == "good",ncol(S_tot_fwd)]))
ddm <- density(log10(S_tot_fwd[gmb == "medium",ncol(S_tot_fwd)]))
ddb <- density(log10(S_tot_fwd[gmb == "bad",ncol(S_tot_fwd)]))
par(mar = c(5.1, 0, 4.1, 0.5))
plot(ddg$y, ddg$x, type = "l", col = cgmb$good, lwd = 2, xlab = "", ylab = "", 
     xlim = c(-0.01,max(ddg$y, ddm$y, ddb$y)*1.05), ylim = par("usr")[3:4], 
     xaxs = "i", xaxt = "n", yaxt = "n")
polygon(c(rep(0,length(ddg$y)), rev(ddg$y)), c(ddg$x, rev(ddg$x)), col = cgmb$good, border = NA)
polygon(c(rep(0,length(ddm$y)), rev(ddm$y)), c(ddm$x, rev(ddm$x)), col = cgmb$medium, border = NA)
polygon(c(rep(0,length(ddb$y)), rev(ddb$y)), c(ddb$x, rev(ddb$x)), col = cgmb$bad, border = NA)

rm(list=c("mod","pop","yy","S_tot","S_tot_fwd","type","log_phi","quantiles","gmb","at",
          "c1","c1t","cgmb","alpha","indx","ddg","ddm","ddb"))
dev.off()


#------------------------------------------------------------------------------------
# Pop-level probabilities of quasi-extinction under multi-pop IPM
# as a function of harvest rate, conditioned on average environmental conditions
# (good/medium/bad)
#------------------------------------------------------------------------------------

dev.new(width = 10, height = 7)
png(filename="PQE_vs_F_given_phi.png", width=10, height=7, units="in", res=200, type="cairo-png")
par(mfrow=c(2,3), mar=c(1,3,4.1,2.1), oma=c(5.1,2.1,2.1,0))

mod <- "PVA_F"
qet <- 50     # set quasi-extinction threshold (4-year moving average)
tfin <- 25    # set time horizon
pops <- c("Catherine","Bear Valley","Secesh","Lemhi","Yankee","Wenatchee")
type <- ifelse(years %in% fish_data$year, "past", "future")
log_phi <- log(eval(as.name(mod))[,,substring(dimnames(eval(as.name(mod)))[[3]],1,4)=="phi["])
log_phi <- matrix(log_phi, nrow = prod(dim(log_phi)[1:2]))
log_phi <- log_phi[,type=="future"][,1:tfin]
quantiles <- quantile(rowMeans(log_phi), c(1/3, 2/3))
gmb <- ifelse(rowMeans(log_phi) >= quantiles[1],
              ifelse(rowMeans(log_phi) >= quantiles[2], "good", "medium"), "bad")
S_tot_fwd <- eval(as.name(mod))[,,substring(dimnames(eval(as.name(mod)))[[3]],1,9)=="S_tot_fwd"]
S_tot_fwd <- matrix(S_tot_fwd, nrow = prod(dim(S_tot_fwd)[1:2]))

for(i in pops)
{
  plot(range(F_rate_fwd), 0:1, pch = "", las = 1, cex.lab = 1.5, cex.axis = 1.2,
       xaxs = "i", yaxs = "i", xlab = "", ylab = "")
  mtext(i, side = 3, line = 0.5, cex = par("cex")*1.5)
  if(i %in% pops[4:6]) mtext("Exploitation rate", side = 1, line = 3.5, cex = par("cex")*2)
  
  pqe_F_g <- vector("numeric", length(F_rate_fwd))
  pqe_F_m <- vector("numeric", length(F_rate_fwd))
  pqe_F_b <- vector("numeric", length(F_rate_fwd))
  for(j in 1:length(F_rate_fwd))
  {
    ss <- S_tot_fwd[gmb == "good",fish_data_fwd$pop==i & fish_data_fwd$F_rate==F_rate_fwd[j]]
    ss <- ss[,1:tfin]
    pqe_F_g[j] <- mean(apply(ss, 1, function(x) any(rollmean(x, 4) < qet)))
    
    ss <- S_tot_fwd[gmb == "medium",fish_data_fwd$pop==i & fish_data_fwd$F_rate==F_rate_fwd[j]]
    ss <- ss[,1:tfin]
    pqe_F_m[j] <- mean(apply(ss, 1, function(x) any(rollmean(x, 4) < qet)))
    
    ss <- S_tot_fwd[gmb == "bad",fish_data_fwd$pop==i & fish_data_fwd$F_rate==F_rate_fwd[j]]
    ss <- ss[,1:tfin]
    pqe_F_b[j] <- mean(apply(ss, 1, function(x) any(rollmean(x, 4) < qet)))
    
    rm(ss)
  }
  
  lines(F_rate_fwd, pqe_F_g, col = "green", lwd = 3)
  lines(F_rate_fwd, pqe_F_m, col = "yellow2", lwd = 3)
  lines(F_rate_fwd, pqe_F_b, col = "red", lwd = 3)
}
mtext(paste("Probability of quasi-extinction in", tfin, "yr"), 
      side = 2, line = 0, outer = T, cex = par("cex")*2)
mtext(ifelse(mod == "PVA_F", "Constant baseline", "Step change 1970"), 
      side = 3, line = 0, outer = T, cex = par("cex")*2)

rm(list = c("qet","pops","pqe_F_g","pqe_F_m","pqe_F_b","type","log_phi","S_tot_fwd",
            "quantiles","gmb","tfin"))
dev.off()


#------------------------------------------------------------------------------------
# Posterior predictive distribution of final population size under multi-pop IPM
# as a function of harvest rate, conditioned on average environmental conditions
# (good/medium/bad)
#------------------------------------------------------------------------------------

dev.new(width = 10, height = 7)
png(filename="S_tot_vs_F_given_phi.png", width=10, height=7, units="in", res=200, type="cairo-png")
par(mfrow=c(2,3), mar=c(1,3,4.1,2.1), oma=c(5.1,2.1,2.1,0))

mod <- "PVA_F"
tfin <- 25    # set time horizon
pops <- c("Catherine","Bear Valley","Secesh","Lemhi","Yankee","Wenatchee")
yfin <- sapply(pops, function(i) max(fish_data$year[fish_data$pop==i]) + tfin)
type <- ifelse(years %in% fish_data$year, "past", "future")
log_phi <- log(eval(as.name(mod))[,,substring(dimnames(eval(as.name(mod)))[[3]],1,4)=="phi["])
log_phi <- matrix(log_phi, nrow = prod(dim(log_phi)[1:2]))
log_phi <- log_phi[,type=="future"][,1:tfin]
quantiles <- quantile(rowMeans(log_phi), c(1/3, 2/3))
gmb <- ifelse(rowMeans(log_phi) >= quantiles[1],
              ifelse(rowMeans(log_phi) >= quantiles[2], "good", "medium"), "bad")
alpha <- c(0.2, 0.3, 0.2)
cgmb <- c(good = "green3", medium = "yellow2", bad = "red2")
cgmbt <- lapply(cgmb, col2rgb)
cgmbt <- sapply(c(good = 1, medium = 2, bad = 3), 
               function(i) rgb(cgmbt[[i]][1], cgmbt[[i]][2], cgmbt[[i]][3], maxColorValue = 255, alpha = 255*alpha[i]))
S_tot_fwd <- eval(as.name(mod))[,,substring(dimnames(eval(as.name(mod)))[[3]],1,9)=="S_tot_fwd"]
S_tot_fwd <- matrix(S_tot_fwd, nrow = prod(dim(S_tot_fwd)[1:2]))

for(i in pops)
{
  FF <- fish_data_fwd$F_rate[fish_data_fwd$pop==i & fish_data_fwd$year==yfin[i]]
  S_tot_F <- S_tot_fwd[,fish_data_fwd$pop==i & fish_data_fwd$year==yfin[i]]

  plot(FF, apply(S_tot_F, 2, median), pch = "", las = 1, cex.lab = 1.5, cex.axis = 1.2,
       ylim = c(5e-2, 1e5),
       # ylim = range(apply(S_tot_fwd[,fish_data_fwd$pop %in% pops & fish_data_fwd$year %in% yfin], 2, quantile, c(0.01,0.99))),
       xaxs = "i", yaxt = "n", log = "y", xlab = "", ylab = "")
  for(j in c("bad","medium","good"))
    polygon(c(FF, rev(FF)), 
             c(apply(S_tot_F[gmb==j,], 2, quantile, 0.025), 
               rev(apply(S_tot_F[gmb==j,], 2, quantile, 0.975))),
             col = cgmbt[j], border = NA)
  for(j in c("bad","medium","good"))
    lines(FF, apply(S_tot_F[gmb==j,], 2, median), lwd = 3, col = cgmb[j])
  at <- maglab(10^par("usr")[3:4], log = T)
  axis(2, at$labat, cex.axis = 1.2, las = 1,
       labels = sapply(log10(at$labat), function(i) as.expression(bquote(10^ .(i)))))
  mtext(i, side = 3, line = 0.5, cex = par("cex")*1.5)
  if(i %in% pops[4:6]) mtext("Exploitation rate", side = 1, line = 3.5, cex = par("cex")*2)
}
mtext(paste("Spawners in year", tfin), side = 2, line = 0, outer = T, cex = par("cex")*2)
mtext(ifelse(mod == "PVA_F", "Constant baseline", "Step change 1970"), 
      side = 3, line = 0, outer = T, cex = par("cex")*2)

rm(list = c("pops","S_tot_F","type","log_phi","S_tot_fwd","cgmb","cgmbt","alpha",
            "quantiles","gmb","tfin","yfin","FF"))
dev.off()


#------------------------------------------------------------------------------------
# Pop-level unconditional probabilities of quasi-extinction under multi-pop IPM
# as a function of harvest rate
#------------------------------------------------------------------------------------

dev.new(width = 7, height = 7)
# png(filename="PQE_vs_F.png", width=7, height=7, units="in", res=200, type="cairo-png")
qet <- 50     # set quasi-extinction threshold (4-yr moving average)
S_tot_fwd <- PVA_F_1970[,,grep("S_tot_fwd", dimnames(PVA_F_1970)[[3]])]
S_tot_fwd <- t(matrix(S_tot_fwd, nrow = prod(dim(S_tot_fwd)[1:2])))
pqe_F <- aggregate(S_tot_fwd, list(pop = fish_data_fwd$pop, F_rate = fish_data_fwd$F_rate), 
                   function(x) any(rollmean(x, 4) < qet))
c2 <- "blue4"
c2t <- col2rgb(c2)
c2t <- rgb(c2t[1], c2t[2], c2t[3], maxColorValue = 255, alpha = 255*0.7)

plot(range(pqe_F$F_rate), 0:1, pch = "", las = 1, cex.lab = 1.5, cex.axis = 1.2, 
     xaxs = "i", yaxs = "i", xlab = "Exploitation rate", ylab = "Probability of quasi-extinction")
for(i in levels(pqe_F$pop))
  lines(pqe_F$F_rate[pqe_F$pop==i], rowMeans(pqe_F[pqe_F$pop==i,-(1:2)]), col = c2t)

rm(list = c("qet","c2","c2t","S_tot_fwd","pqe_F"))
# dev.off()





