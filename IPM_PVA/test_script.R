# load data
fish_data <- read.table(file.path("~", "SalmonIPM", "IPM_PVA", "fish_data.txt"), sep = "\t", header = T)

# fit hierarchical spawner-recruit model to run reconstruction
test_fit <- salmonIPM(fish_data = fish_data, model = "RR", chains = 3, iter = 1000, warmup = 500)

print(test_fit, pars = c("phi","R_hat"), include = F)
