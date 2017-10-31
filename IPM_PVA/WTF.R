windows()
par(mfrow=c(3,1))

dlp <- extract(IPM_pp2b,"lp__", permute = F, inc_warmup = T) - 
  extract(IPM_pp2,"lp__", permute = F, inc_warmup = T)
for(i in 1:3)
{
  plot(1:1000, dlp[,i,1], pch="", xaxs = "i", 
       xlab = "iter", ylab = "difference in lp__", main = paste0("chain ", i)) 
  abline(v = 500, col = "gray")
  abline(h = 0, lty = 2)
  lines(1:1000, dlp[,i,1], col = "blue")
}
rm(dlp)


windows()
plot(summary(IPM_pp2b)$summary[,"n_eff"], summary(IPM_pp2)$summary[,"n_eff"],
     xlim = c(0,1500), ylim = c(0,1500), xlab = "model 1", ylab = "model 2", main = "n_eff")
abline(0,1)

neff1 <- summary(IPM_pp2b)$summary[,"n_eff"]
neff2 <-  summary(IPM_pp2)$summary[,"n_eff"]
dneff <- neff1 - neff2
summary(dneff[!(neff1 == 1500 & neff2 == 1500)])
t.test(dneff[!(neff1 == 1500 & neff2 == 1500)])
rm(dneff);rm(neff1);rm(neff2)