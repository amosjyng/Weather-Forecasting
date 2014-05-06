library(glasso) # keep an eye out for dpglasso
library(hydroGOF)

# read data from commandline argument
args <- commandArgs(trailingOnly = TRUE)
d <- read.delim(args[1], header = FALSE)
#colnames(d) <- c("P", "PE", "CP", "LW", "SW", "PE", "SP", "SH", "T", "WS")
colnames(d) <- c("P", "PE", "CP", "LW", "SW", "PE", "SP", "SH", "T", "WS", "P2", "PE2", "CP2", "LW2", "SW2", "PE2", "SP2", "SH2", "T2", "WS2")
#colnames(d) <- c("P", "PE", "CP", "LW", "SW", "PE", "SP", "SH", "T", "WS", "gP", "gPE", "gCP", "gLW", "gSW", "gPE", "gSP", "gSH", "gT", "gWS", "P2", "PE2", "CP2", "LW2", "SW2", "PE2", "SP2", "SH2", "T2", "WS2")
regularization <- as.numeric(args[2])

predict <- function(i, gm, train, test) {
  mus <- colMeans(d)
  vars <- gm$w
  # condition  on everything but the variable in the i-th column
  nvars <- ncol(train)
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
  
  # test
  muBar <- numeric(nrow(test))
  sigmaBar <- sigma11 - sigma12 %*% solve(sigma22) %*% sigma21
  for (row in 1:nrow(test)) {
    muBar[[row]] <- mu1 + sigma12 %*% solve(sigma22) %*% (t(as.matrix(test[row,allExcept_i])) - mu2)
  }
  rmse <- rmse(muBar, test[,i])
  
  return(c(rmse, sigmaBar))
}

test <- function(train, test) {
  stats <- matrix(-1, nrow = 2, ncol = ncol(train))
  print("Training...")
  gm <- glasso(var(train), regularization)
  
  print("Testing...")
  for (i in 1:ncol(train)) {
    print(cat(sprintf("Predicting variable #%i\n", i)))
    results <- predict(i, gm, train, test)
    stats[1,i] <- results[1] / sd(test[,i])
    stats[2,i] <- results[2] / gm$w[i,i]
  }
  
  rownames(stats) <- c("RMSE / St. Dev", "Sigma Bar / Sigma")
  #colnames(stats) <- c("P", "PE", "CP", "LW", "SW", "PE", "SP", "SH", "T", "WS")
  colnames(stats) <- c("P", "PE", "CP", "LW", "SW", "PE", "SP", "SH", "T", "WS", "P2", "PE2", "CP2", "LW2", "SW2", "PE2", "SP2", "SH2", "T2", "WS2")
  #colnames(stats) <- c("P", "PE", "CP", "LW", "SW", "PE", "SP", "SH", "T", "WS", "gP", "gPE", "gCP", "gLW", "gSW", "gPE", "gSP", "gSH", "gT", "gWS", "P2", "PE2", "CP2", "LW2", "SW2", "PE2", "SP2", "SH2", "T2", "WS2")
  
  return(stats)
}

total <- nrow(d)
chunk <- as.integer(total / 10)

cross_validate <- function(i) {
  start <- (i - 1) * chunk + 1
  end <- i * chunk
  upto <- NULL
  if (start > 1) {
    upto <- 1:(start - 1)
  }
  from <- NULL
  if (end < nrow(d)) {
    from <- (end + 1):nrow(d)
  }
  test <- d[start:end,]
  train <- d[c(upto,from),]
  
  stats <- test(train, test)
  print(cat(sprintf("For cross-validation iteration %i:", i)))
  print(stats)
  
  return (stats)
}

stat_matrices <- list()
for (i in 1:10) {
  stat_matrices[[i]] <- cross_validate(i)
}

add <- function(x) Reduce("+", x)

avg_stats <- add(stat_matrices) / 10
print(avg_stats)
save(avg_stats, file="stats.Rda")