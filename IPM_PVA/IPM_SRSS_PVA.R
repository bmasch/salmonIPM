setwd(file.path("~","salmonIPM","IPM_PVA"))
options(device=windows)
library(salmonIPM)
library(magicaxis)
library(zoo)

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
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.25)
c1tt <- rgb(c1tt[1], c1tt[2], c1tt[3], maxColorValue = 255, alpha = 255*0.45)
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
          col = c1t, border = NA)
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
          col = c1tt, border = NA)
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
BH <- function(a, b, S) 
{
  a*S/(1 + b*S)
}

S <- matrix(seq(0, quantile(fish_data$S_tot_obs/fish_data$A, 0.9, na.rm = T), length = 100),
            nrow = sum(PVA_RR_pp@sim$n_save - PVA_RR_pp@sim$warmup2), ncol = 100, byrow = T)

# S-R curves
mu_log_a <- as.vector(extract1(PVA_RR_pp,"mu_log_a"))
mu_log_b <- as.vector(extract1(PVA_RR_pp,"mu_log_b"))
R_ESU_RR <- BH(a = exp(mu_log_a), b = exp(mu_log_b), S = S)
mu_log_a <- as.vector(extract1(PVA_IPM_pp,"mu_log_a"))
mu_log_b <- as.vector(extract1(PVA_IPM_pp,"mu_log_b"))
R_ESU_IPM <- BH(a = exp(mu_log_a), b = exp(mu_log_b), S = S)

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
dd_IPM_ESU <- density(extract1(PVA_IPM_pp,"mu_log_a") - extract1(PVA_IPM_pp,"mu_log_b"))
dd_IPM_pop <- vector("list", length(levels(fish_data$pop)))
for(i in 1:length(dd_IPM_pop))
  dd_IPM_pop[[i]] <- density(log(extract1(PVA_IPM_pp,"a")[,i]) - log(extract1(PVA_IPM_pp,"b")[,i]))
dd_RR_ESU <- density(extract1(PVA_RR_pp,"mu_log_a") - extract1(PVA_RR_pp,"mu_log_b"))
dd_RR_pop <- vector("list", length(levels(fish_data$pop)))
for(i in 1:length(dd_RR_pop))
  dd_RR_pop[[i]] <- density(log(extract1(PVA_RR_pp,"a")[,i]) - log(extract1(PVA_RR_pp,"b")[,i]))

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

rm(list=c("mu_log_a","mu_log_b","S","R_ESU_RR","R_ESU_IPM","BH",
          "bb","dd_IPM_ESU","dd_RR_ESU","dd_IPM_pop","dd_RR_pop"))
# dev.off()


#-------------------------------------------------------------------------
# Pop-specific S-R curves with and without obs error
#-------------------------------------------------------------------------

dev.new(width = 7, height = 7)
# png(filename="Fig_2.5(3).png", width=7, height=7, units="in", res=200, type="cairo-png")
BH <- function(a, b, S, A) 
{
  a*S/(A + b*S)
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
b <- as.vector(extract1(PVA_RR_pp,"b")[,which(levels(fish_data$pop)==pop)])
R_RR <- BH(a = a, b = b, S = S, AA) * AA
a <- as.vector(extract1(PVA_IPM_pp,"a")[,which(levels(fish_data$pop)==pop)])
b <- as.vector(extract1(PVA_IPM_pp,"b")[,which(levels(fish_data$pop)==pop)])
R_IPM <- BH(a = a, b = b, S = S, AA) * AA

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

rm(list = c("yy","AA","BH","pop","S","a","b","R_RR","R_IPM","S_tot_RR","S_tot_IPM",
            "R_adj_RR","R_adj_IPM","c1","c1t","c2","c2t","c2tt"))

# dev.off()


#------------------------------------------------------------------
# Probability of quasi-extinction by population under RR and IPM
#------------------------------------------------------------------

dev.new(height = 7, width = 7)
# png(filename="Fig_3.png", width=7, height=7, units="in", res=200, type="cairo-png")
par(oma = c(0,5,0,0))
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
     xlab = "Harvest rate", ylab = "Probability of decline")
lines(sort(Umax_ESU_IPM), (1:M)/M, lwd = 4, col = c2)
for(i in 1:ncol(Umax_pop_IPM))
{
  lines(sort(Umax_pop_RR[,i]), (1:M)/M, lwd = 1, col = c1t)
  lines(sort(Umax_pop_IPM[,i]), (1:M)/M, lwd = 1, col = c2t)
}
legend("topleft", c("IPM","RR"), col = c("blue4","orangered3"), lwd = 3, cex = 1.2)

par(mar = c(0.1,3.1,0.5,1.5))
bins <- seq(0, 1, 0.05)
F1 <- hist(fish_data$F_rate[fish_data$year < 1980], breaks = bins, plot = F)$counts
F2 <- hist(fish_data$F_rate[fish_data$year >= 1980], breaks = bins, plot = F)$counts
barplot(rbind(F1,F2), names = bins[-1], space = 0, yaxs = "i", xaxt = "n", yaxt = "n", 
     xlab = "", ylab = "", main = "", col = c("black","gray"), border = "white")

rm(list = c("mu_log_a_RR","mu_log_a_IPM","Umax_ESU_RR","Umax_ESU_IPM",
            "a_RR","Umax_pop_RR","a_IPM","Umax_pop_IPM","c1","c1t","c2","c2t","bins","F1","F2"))

# dev.off()









#===========================================================================
# MORE FIGURES: SPARE PARTS
#===========================================================================

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
# Time series of observed and fitted or predicted total spawners for each pop
#--------------------------------------------------------------------------------

dev.new(width=16,height=10)
par(mfrow=c(4,6), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))
S_tot <- extract1(PVA_IPM_pp,"S_tot")
c1 <- "darkgray"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.6)
c2 <- "blue"
c2t <- col2rgb(c2)
c2t <- rgb(c2t[1], c2t[2], c2t[3], maxColorValue = 255, alpha = 255*0.4)
iters <- 1:5
for(i in levels(fish_data_aug$pop))
{
  S_tot_i <- S_tot[,fish_data_aug$pop==i]
  S_tot_fit_i <- S_tot_i[,fish_data_aug$type[fish_data_aug$pop==i]=="past"]
  # S_tot_fit_pts_i <- density(log(S_tot_fit_i[,1]), 
  #                          from = min(apply(log(S_tot_i), 2, quantile, 0.025)), 
  #                          to = max(apply(log(S_tot_i), 2, quantile, 0.975)))$x
  # S_tot_fit_dens_i <- apply(log(S_tot_fit_i), 2, function(x) 
  #   density(x, from = min(apply(log(S_tot_i), 2, quantile, 0.025)), 
  #           to = max(apply(log(S_tot_i), 2, quantile, 0.975)))$y)
  S_tot_fore_i <- cbind(S_tot_fit_i[,ncol(S_tot_fit_i)],
                        S_tot_i[,fish_data_aug$type[fish_data_aug$pop==i]=="future"])
  # S_tot_fore_pts_i <- density(log(S_tot_fore_i[,1]), 
  #                            from = min(apply(log(S_tot_i), 2, quantile, 0.025)), 
  #                            to = max(apply(log(S_tot_i), 2, quantile, 0.975)))$x
  # S_tot_fore_dens_i <- apply(log(S_tot_fore_i), 2, function(x) 
  #   density(x, from = min(apply(log(S_tot_i), 2, quantile, 0.025)), 
  #           to = max(apply(log(S_tot_i), 2, quantile, 0.975)))$y)
  year_fit_i <- fish_data_aug$year[fish_data_aug$pop==i & fish_data_aug$type=="past"]
  year_fore_i <- c(year_fit_i[length(year_fit_i)], 
                   fish_data_aug$year[fish_data_aug$pop==i & fish_data_aug$type=="future"])
  
  plot(fish_data_aug$year[fish_data_aug$pop==i], fish_data_aug$S_tot_obs[fish_data_aug$pop==i], pch="", cex=1.2, cex.axis=1.2, las=1,
       # ylim = c(0,max(apply(S_tot_i, 2, quantile, 0.975))), xlab="", ylab="")
       # ylim = c(min(apply(S_tot_i, 2, quantile, 0.025)),max(apply(S_tot_i, 2, quantile, 0.975))),
       xlim = range(fish_data_aug$year),
       ylim = range(c(range(S_tot_i[iters,]),
                      c(min(apply(S_tot_fit_i, 2, quantile, 0.025)),max(apply(S_tot_fit_i, 2, quantile, 0.975))),
                      range(replace(fish_data_aug$S_tot_obs[fish_data_aug$pop==i], fish_data_aug$S_tot_obs[fish_data_aug$pop==i]==0, 1), 
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
  points(fish_data_aug$year[fish_data_aug$pop==i], fish_data_aug$S_tot_obs[fish_data_aug$pop==i], type="b", pch=16)
  mtext(i, side=3, line=0.5, cex=1.2)
}
mtext("Year", outer = T, side=1, line=2, cex=1.2)
mtext("Spawners", outer = T, side=2, line=1.1, cex=1.2)
rm(list=c("mod","S_tot","S_tot_i","S_tot_fit_i","S_tot_fore_i","year_fit_i","year_fore_i",
          "c1","c1t","c2","c2t","at","iters"))


#-----------------------------------------------------------------------------------------
# Time series of observed and fitted or predicted recruits per spawner for each pop
#-----------------------------------------------------------------------------------------

dev.new(width=16,height=10)
par(mfrow=c(4,6), mar=c(1,2,4.1,1), oma=c(4.1,3.1,0,0))
S_tot <- extract1(PVA_IPM_pp,"S_tot")
R_tot <- extract1(PVA_IPM_pp,"R_tot")
RS <- R_tot/S_tot
rr <- run_recon(fish_data_aug)
c1 <- "darkgray"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.6)
c2 <- "blue"
c2t <- col2rgb(c2)
c2t <- rgb(c2t[1], c2t[2], c2t[3], maxColorValue = 255, alpha = 255*0.4)
iters <- 1:5
for(i in levels(fish_data_aug$pop))
{
  RS_obs_i <- (rr$R/rr$S)[rr$pop==i]
  RS_i <- RS[,fish_data_aug$pop==i]
  RS_fit_i <- RS_i[,fish_data_aug$type[fish_data_aug$pop==i]=="past"]
  RS_fore_i <- cbind(RS_fit_i[,ncol(RS_fit_i)],
                     RS_i[,fish_data_aug$type[fish_data_aug$pop==i]=="future"])
  year_fit_i <- fish_data_aug$year[fish_data_aug$pop==i & fish_data_aug$type=="past"]
  year_fore_i <- c(year_fit_i[length(year_fit_i)], 
                   fish_data_aug$year[fish_data_aug$pop==i & fish_data_aug$type=="future"])
  
  plot(fish_data_aug$year[fish_data_aug$pop==i], RS_obs_i, pch = "", cex=1.2, cex.axis=1.2, las=1,
       # ylim=c(min(apply(RS_i, 2, quantile, 0.025)),max(apply(RS_i, 2, quantile, 0.975))), 
       xlim = range(fish_data_aug$year),
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
  points(fish_data_aug$year[fish_data_aug$pop==i], RS_obs_i, type = "b", pch = 16)
  mtext(i, side=3, line=0.5, cex=1.2)
}
mtext("Year", outer = T, side=1, line=2, cex=1.2)
mtext("Recruits / spawner", outer = T, side=2, line=1.1, cex=1.2)
rm(list=c("rr","S_tot","RS_obs_i","RS","RS_i","RS_fit_i","RS_fore_i",
          "year_fit_i","year_fore_i","c1","c1t","c2","c2t","at","iters"))


#--------------------------------------------------
# Shared brood-year productivity anomalies
#--------------------------------------------------

dev.new(height = 10, width = 10)
# png(filename="IPM_SRSS_PVA_year_effects.png", width=10, height=10, units="in", res=200, type="cairo-png")
par(mar = c(5.1,5.1,4.1,2.1))
phi <- extract1(PVA_IPM_pp, "phi")
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











