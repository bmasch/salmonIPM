# load data
fish_data <- read.table(path.expand(file.path("~", "SalmonIPM", "tmp_data", "fish_data.txt")), sep = "\t", header = T)

# fit hierarchical run reconstruction spawner-recruit model
test_fit <- salmonIPM(fish_dat = fish_data, model = "RR", chains = 3, iter = 1000, warmup = 500)

print(test_fit, pars = c("phi","R_hat"), include = FALSE)
