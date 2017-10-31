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
plot(summary(IPM_pp1)$summary[,"n_eff"], summary(IPM_pp2)$summary[,"n_eff"])
abline(0,1)
