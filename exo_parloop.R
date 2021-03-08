set.seed(020193)

crs <- detectCores()
cl <- makeCluster(getOption("cl.cores", crs-1))
writeLines(paste("Running with", crs-1, 'cores'))
# Load packages to cluster
invisible(clusterCall(cl, function(pkgs) {
  # library(MASS)
  # library(forecast)
  # library(tsutils)
  library(clusterGeneration)
  library(smooth)
  # library(nloptr)
}))
invisible(clusterExport(cl, "iter.run"))
invisible(clusterExport(cl, "myLoss_BackcastingInitials"))
invisible(clusterExport(cl, "myLoss_OptimalInitials"))

runs <- 500

system.time({resultEXO <- clusterApplyLB(cl, 1:runs, function(x) iter.run())})

stopCluster(cl)


set.seed(020193)
mat.benchmark <- matrix(NA, ncol = 2, nrow = 500)
for (i in 1:500) {
  mat.benchmark[i,] <- iter.run.adam()
}
