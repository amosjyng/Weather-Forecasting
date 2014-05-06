library(glasso) # keep an eye out for dpglasso
library(hydroGOF)

# read data from commandline argument
args <- commandArgs(trailingOnly = TRUE)
data <- read.delim(args[1], header = TRUE)

# run GLASSO algorithm
load("glasso_model.Rda")
mus <- colMeans(data)
vars <- gm$w

# condition  on everything but p
mu1 <- as.matrix(mus[1])
mu2 <- as.matrix(mus[2:10])
sigma11 <- as.matrix(vars[1,1])
sigma12 <- t(as.matrix(vars[1,2:10]))
sigma21 <- as.matrix(vars[2:10,1])
sigma22 <- as.matrix(vars[2:10,2:10])

muBar <- numeric(nrow(data))
sigmaBar <- sigma11 - sigma12 %*% solve(sigma22) %*% sigma21
for (row in 1:nrow(data)) {
  muBar[[row]] <- mu1 + sigma12 %*% solve(sigma22) %*% (t(as.matrix(data[row,2:10])) - mu2)
}
rmse <- rmse(muBar, data[,1])
stdev <- sd(data[,1])

print(cat(sprintf("RMSE / St. Dev = %f\n", rmse / stdev)))
print(cat(sprintf("Sigma Bar / Sigma = %f\n", sigmaBar / vars[1,1])))