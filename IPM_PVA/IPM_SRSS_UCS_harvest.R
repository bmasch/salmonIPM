setwd(file.path("~","salmonIPM","IPM_PVA"))
options(device=windows)
library(salmonIPM)

#===========================================================================
# DATA
#===========================================================================

# Load data
fish_data <- read.table(file.path("~", "salmonIPM", "IPM_PVA", "fish_data.txt"), sep = "\t", header = T)
fish_data <- fish_data[order(fish_data$code, fish_data$year),]

# Pad data with years through max_year
N_future_years <- 50
max_year <- max(fish_data$year) + N_future_years
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
# HARVEST SIMULATIONS
# Vary future exploitation rate (fixed-rate policy) and look at
# impact on quasi-extinction risk under multi-pop IPM
#===========================================================================

n_chains <- 3
n_warmup <- 500
n_save <- 500

F_rate_future <- seq(0, 0.3, length=7)
S_tot_F <- array(NA, dim = c(n_chains*n_save, nrow(fish_data_aug), length(F_rate_future)))
R_tot_F <- array(NA, dim = c(n_chains*n_save, nrow(fish_data_aug), length(F_rate_future)))
log_phi_F <- array(NA, dim = c(n_chains*n_save, length(unique(fish_data_aug$year)), length(F_rate_future)))

for(i in 1:length(F_rate_future))
{
  fish_data_F <- fish_data_aug
  fish_data_F$F_rate[fish_data_F$type=="future"] <- F_rate_future[i]
  
  PVA_F <- salmonIPM(fish_data = fish_data_F, model = "IPM", pool_pops = TRUE, 
                     pars = c("R_tot","S_tot"),
                     chains = n_chains, iter = n_warmup + n_save, warmup = n_warmup,
                     control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))
  S_tot_F[,,i] <- extract1(PVA_F,"S_tot")
  R_tot_F[,,i] <- extract1(PVA_F,"R_tot")
  log_phi_F[,,i] <- log(extract1(PVA_F,"phi"))
}

rm(list = c("n_chains","n_warmup","n_save"))


#===========================================================================
# FIGURES
#===========================================================================

#------------------------------------------------------------------------------------
# Pop-level unconditional probabilities of quasi-extinction under multi-pop IPM
# as a function of harvest rate
#------------------------------------------------------------------------------------

dev.new(width = 7, height = 7)
# png(filename="PQE_vs_F.png", width=7, height=7, units="in", res=200, type="cairo-png")
qet <- 50     # set quasi-extinction threshold (4-yr moving average)
pop <- fish_data_F$pop[fish_data_F$type=="future"]
c2 <- "blue4"
c2t <- col2rgb(c2)
c2t <- rgb(c2t[1], c2t[2], c2t[3], maxColorValue = 255, alpha = 255*0.7)

plot(range(F_rate_future), 0:1, pch = "", las = 1, cex.lab = 1.5, cex.axis = 1.2, 
     xaxs = "i", yaxs = "i", xlab = "Exploitation rate", ylab = "Probability of quasi-extinction")
for(i in levels(pop))
{
  pqe_F <- S_tot_F[,fish_data_F$pop==i & fish_data_F$type=="future",]
  pqe_F <- colMeans(apply(pqe_F, c(1,3), function(x) any(rollmean(x, 4) < qet)))
  lines(F_rate_future, pqe_F, col = c2t)
}

rm(list = c("qet","pop","c2","c2t","pqe_F"))
# dev.off()


#------------------------------------------------------------------------------------
# Sample paths of brood-year log productivity anomalies (log phi)
#------------------------------------------------------------------------------------

dev.new(width = 14, height = 7)
# png(filename="log_phi_sample_paths.png", width=14, height=7, units="in", res=200, type="cairo-png")
layout(matrix(c(1,2), nrow = 1), widths = c(14,2))

yr <- sort(unique(fish_data_aug$year))
type <- ifelse(yr %in% fish_data$year, "past", "future")
log_phi <- log_phi_F[,,1]
quantiles <- quantile(rowMeans(log_phi[,type == "future"]), c(1/3, 2/3))

c1 <- "blue4"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.5)

indx1 <- c(which.max(rowMeans(log_phi[,type=="future"])), which.min(rowMeans(log_phi[,type=="future"])))
indx2 <- sample(nrow(log_phi),100)

# layer 1
par(mar = c(5.1, 4.5, 4.1, 0.2))
plot(yr, log_phi[indx1[1],], type = "l", col = c1t,
     las = 1, cex.lab = 1.8, cex.axis = 1.5, xaxs = "i", yaxs = "i", ylim = range(log_phi[indx1,]), 
     xlab = "Brood year", ylab = "Productivity anomaly")
abline(v = max(yr[type == "past"]), lty = 2)

# layer 2
lines(yr, log_phi[indx1[2],], col = c1t)

# layer 3
segments(x0 = max(yr[type == "past"]), x1 = max(yr), y0 = rowMeans(log_phi[indx1,]), col = c1t)

# layer 4
for(i in indx2)
  lines(yr, log_phi[i,], col = c1t)
dd <- density(rowMeans(log_phi))
par(mar = c(5.1, 0, 4.1, 0.5))
plot(dd$y, dd$x, type = "l", ylim = par("usr")[3:4], col = c1, lwd = 2,
     xlab = "", ylab = "", xlim = c(-0.1,max(dd$y)*1.05), xaxs = "i", xaxt = "n", yaxt = "n")
x1 <- which.min(abs(dd$x - quantiles[1]))
x2 <- which.min(abs(dd$x - quantiles[2]))
x3 <- length(dd$x)
polygon(c(rep(0,x1), dd$y[1:x1]), c(dd$x[x1:1], dd$x[1:x1]), col = "red", border = NA)
polygon(c(rep(0, x2 - x1), dd$y[(x1+1):x2]), c(dd$x[x2:(x1+1)], dd$x[(x1+1):x2]), col = "yellow", border = NA)
polygon(c(rep(0, x3 - x2), dd$y[(x2+1):x3]), c(dd$x[x3:(x2+1)], dd$x[(x2+1):x3]), col = "green", border = NA)

segments(x0 = 0, x1 = dd$y[which.min(abs(dd$x - quantiles[1]))], y0 = quantiles[1], col = c1)
segments(x0 = 0, x1 = dd$y[which.min(abs(dd$x - quantiles[2]))], y0 = quantiles[2], col = c1)

rm(list = c("yr","type","log_phi","c1","c1t","indx1","indx2","dd","quantiles","x1","x2","x3"))
# dev.off()


#------------------------------------------------------------------------------------
# Pop-level probabilities of quasi-extinction under multi-pop IPM
# as a function of harvest rate, conditioned on average environmental conditions
# (good/medium/bad)
#------------------------------------------------------------------------------------

# dev.new(width = 7, height = 7)
png(filename="PQE_vs_F_given_phi.png", width=7, height=7, units="in", res=200, type="cairo-png")
par(mfrow=c(2,2), mar=c(1,3,4.1,2.1), oma=c(5.1,2.1,0,0))

qet <- 50     # set quasi-extinction threshold (4-yr moving average)
pops <- c("Bear Valley","Lemhi","Yankee","Wenatchee")

yr <- sort(unique(fish_data_aug$year))
type <- ifelse(yr %in% fish_data$year, "past", "future")
gmb <- array(0, dim = dim(log_phi_F)[c(1,3)])
for(j in 1:length(F_rate_future))
{
  quantiles <- quantile(rowMeans(log_phi_F[,type == "future",j]), c(1/3, 2/3))
  gmb[,j] <- ifelse(rowMeans(log_phi_F[,type == "future",j]) >= quantiles[1],
                ifelse(rowMeans(log_phi_F[,type == "future",j]) >= quantiles[2], "good", "medium"),
                "bad")
}

for(i in pops)
{
  plot(range(F_rate_future), 0:1, pch = "", las = 1, cex.lab = 1.5, cex.axis = 1.2, 
       xaxs = "i", yaxs = "i", main = i, xlab = "", ylab = "")
  if(i %in% pops[c(3,4)]) mtext("Exploitation rate", side = 1, line = 3.5, cex = par("cex")*2)
  
  pqe_F_g <- vector("numeric", length(F_rate_future))
  pqe_F_m <- vector("numeric", length(F_rate_future))
  pqe_F_b <- vector("numeric", length(F_rate_future))
  S_tot_F_i <- S_tot_F[,fish_data_F$pop==i & fish_data_F$type=="future",]
  for(j in 1:length(F_rate_future))
  {
    ss <- S_tot_F[gmb[,j] == "good",fish_data_F$pop==i & fish_data_F$type=="future",j]
    pqe_F_g[j] <- mean(apply(ss, 1, function(x) any(rollmean(x, 4) < qet)))
    
    ss <- S_tot_F[gmb[,j] == "medium",fish_data_F$pop==i & fish_data_F$type=="future",j]
    pqe_F_m[j] <- mean(apply(ss, 1, function(x) any(rollmean(x, 4) < qet)))
    
    ss <- S_tot_F[gmb[,j] == "bad",fish_data_F$pop==i & fish_data_F$type=="future",j]
    pqe_F_b[j] <- mean(apply(ss, 1, function(x) any(rollmean(x, 4) < qet)))
    
    rm(ss)
  }
  
  lines(F_rate_future, pqe_F_g, col = "green", lwd = 3)
  lines(F_rate_future, pqe_F_m, col = "yellow2", lwd = 3)
  lines(F_rate_future, pqe_F_b, col = "red", lwd = 3)
}
mtext("Probability of quasi-extinction", side = 2, line = 0, outer = T, cex = par("cex")*2)

rm(list = c("qet","pops","pqe_F_g","pqe_F_m","pqe_F_b","yr","type","log_phi",
            "quantiles","gmb","S_tot_F_i"))
dev.off()





