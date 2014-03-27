library(glasso) # keep an eye out for dpglasso

# read data from commandline argument
args <- commandArgs(trailingOnly = TRUE)
data <- read.delim(args[1])
s <- var(data)

# run GLASSO algorithm
gm <- glasso(s, 50)

# do something with results
save(gm, file="glasso_model.Rda")
print(paste("Log-likelihoode of model is", toString(gm$loglik)))