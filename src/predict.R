library(glasso) # keep an eye out for dpglasso
library(hydroGOF)

# read data from commandline argument
args <- commandArgs(trailingOnly = TRUE)
data <- read.delim(args[1], header = TRUE)

# load saved GLASSO model
load("glasso_model.Rda")
mus <- colMeans(data)
vars <- gm$w

predict <- function(i) {
  # condition  on everything but the variable in the i-th column
  nvars <- ncol(data)
  upto <- NULL
  if (i > 1) {
    upto <- 1:(i - 1)
  }
  from <- NULL
  if(i < nvars) {
    from <- (i + 1):nvars
  }
  allExcept_i <- c(upto, from)
  
  mu1 <- as.matrix(mus[i])
  mu2 <- as.matrix(mus[allExcept_i])
  sigma11 <- as.matrix(vars[i,i])
  sigma12 <- t(as.matrix(vars[i,allExcept_i]))
  sigma21 <- as.matrix(vars[allExcept_i,i])
  sigma22 <- as.matrix(vars[allExcept_i,allExcept_i])
  
  muBar <- numeric(nrow(data))
  sigmaBar <- sigma11 - sigma12 %*% solve(sigma22) %*% sigma21
  for (row in 1:nrow(data)) {
    muBar[[row]] <- mu1 + sigma12 %*% solve(sigma22) %*% (t(as.matrix(data[row,allExcept_i])) - mu2)
  }
  rmse <- rmse(muBar, data[,i])
  
  return(c(rmse, sigmaBar))
}

stats <- matrix(-1, nrow = 2, ncol = ncol(data))

for (i in 1:ncol(data)) {
  print(cat(sprintf("Predicting variable #%i\n", i)))
  results <- predict(i)
  stats[1,i] <- results[1] / sd(data[,i])
  stats[2,i] <- results[2] / vars[i,i]
}

rownames(stats) <- c("RMSE / St. Dev", "Sigma Bar / Sigma")
colnames(stats) <- colnames(gm$w)
print(stats)