library(glasso)

# read data from commandline argument
L1_penalty <- 1
args <- commandArgs(trailingOnly = TRUE)
data <- read.delim(args[1])
s <- var(data)

# run GLASSO algorithm
gm <- glasso(s, L1_penalty)

# do something with results
save(gm, file="glasso_model.Rda")
print(gm$loglik)