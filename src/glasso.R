library(glasso) # keep an eye out for dpglasso
library(igraph)

# read data from commandline argument
args <- commandArgs(trailingOnly = TRUE)
data <- read.delim(args[1], header = TRUE)

# run GLASSO algorithm
gm <- glasso(var(data), as.numeric(args[2]))
#colnames(gm$w) <- c("P", "PE", "CP", "LW", "SW", "PE", "SP", "SH", "T", "WS")
colnames(gm$wi) <- c("P", "PE", "CP", "LW", "SW", "PE", "SP", "SH", "T", "WS", "P2", "PE2", "CP2", "LW2", "SW2", "PE2", "SP2", "SH2", "T2", "WS2")
# colnames(gm$wi) <- c("P", "PE", "CP", "LW", "SW", "PE", "SP", "SH", "T", "WS",
#                      "gP", "gPE", "gCP", "gLW", "gSW", "gPE", "gSP", "gSH", "gT",
#                      "gWS", "P2", "PE2", "CP2", "LW2", "SW2", "PE2", "SP2",
#                      "SH2", "T2", "WS2")
rownames(gm$w) <- colnames(gm$w)
colnames(gm$w) <- colnames(gm$w)

# do something with results
save(gm, file="glasso_model.Rda")
print(paste("Log-likelihoode of model is", toString(gm$loglik)))
# plot adjacency matrix
adj <- gm$w
for (row in 1:nrow(adj)) {
  for (col in 1:ncol(adj)) {
    if (adj[row,col] != 0) {
      adj[row,col] = 1
    }
    if (row == col) {
      adj[row, col] = 0
    }
  }
}
png("plot.png")
plot.igraph(as.undirected(graph.adjacency(adj)), vertex.size = 20,
            vertex.label.cex = 1)
dev.off()