set.seed(020193)

crs <- detectCores()
cl <- makeCluster(getOption("cl.cores", crs-1))
writeLines(paste("Running with", crs-1, 'cores'))
# Load packages to cluster
invisible(clusterCall(cl, function(pkgs) {
  library(MASS)
  library(forecast)
  library(tsutils)
  library(clusterGeneration)
  library(smooth)
  library(nloptr)
}))
invisible(clusterExport(cl, "matrixMaker"))
invisible(clusterExport(cl, "gen.forecast"))
invisible(clusterExport(cl, "createInit"))
invisible(clusterExport(cl, "lossFun"))
invisible(clusterExport(cl, "loopLambdaCV"))
invisible(clusterExport(cl, "fun.iter"))
invisible(clusterExport(cl, "freq.ts"))
invisible(clusterExport(cl, "benchmark.ets"))

sqq <- c(0, rev(exp(seq.int(from = log(1), to = log(1e-04), by = -(log(1)-log(1e-04))/(30-1)))))
invisible(clusterExport(cl, "sqq"))

runs <- 10

system.time({res.exp1.trial <- clusterApplyLB(cl, 1:runs, function(x) fun.iter(sq = sqq))})

stopCluster(cl)

save(res.exp1, file = "res.exp1.RData")

res.exp1[[1]][[1]][1,]

res.sc1.outRMSE <- matrix(NA, ncol = 30, nrow = 500)
for (i in 1:30) {
  res.sc1.outRMSE[,i] <- t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[1]]))[,8]
}

res.sc2.outRMSE <- matrix(NA, ncol = 30, nrow = 500)
for (i in 1:30) {
  res.sc2.outRMSE[,i] <- t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[2]]))[,6]
}

res.sc3.outRMSE <- matrix(NA, ncol = 30, nrow = 500)
for (i in 1:30) {
  res.sc3.outRMSE[,i] <- t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[3]]))[,8]
}

res.sc4.outRMSE <- matrix(NA, ncol = 30, nrow = 500)
for (i in 1:30) {
  res.sc4.outRMSE[,i] <- t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[4]]))[,6]
}

boxplot(res.sc1.outRMSE, main = "RMSE in Validation Set \n LASSO with optimal initial values",
        ylab = "Out-RMSE", xlab = "Lambda", xaxt = "n")
lines(colMeans(res.sc1.outRMSE), col = "red", type = "o", pch = 19, bg = "red")
axis(1,at=1:30,labels=round(sqq, 4),las=2, cex.axis = 0.8)

boxplot(res.sc2.outRMSE, main = "RMSE in Validation Set \n LASSO with backcast initial values",
        ylab = "Out-RMSE", xlab = "Lambda", xaxt = "n", ylim = c(0,30))
lines(colMeans(res.sc2.outRMSE), col = "red", type = "o", pch = 19, bg = "red")
axis(1,at=1:30,labels=round(sqq, 4),las=2, cex.axis = 0.8)

boxplot(res.sc3.outRMSE, main = "RMSE in Validation Set \n Adaptive with optimal initial values",
        ylab = "Out-RMSE", xlab = "Lambda", xaxt = "n")
lines(colMeans(res.sc3.outRMSE), col = "red", type = "o", pch = 19, bg = "red")
axis(1,at=1:30,labels=round(sqq, 4),las=2, cex.axis = 0.8)

boxplot(res.sc4.outRMSE, main = "OutRMSE - Adaptive with backcast initial values",
        ylab = "Out-RMSE", xlab = "Lambda", xaxt = "n")
lines(colMeans(res.sc4.outRMSE), col = "red", type = "o", pch = 19, bg = "red")
axis(1,at=1:30,labels=round(sqq,4),las=2, cex.axis = 0.8)


res.sc1.smoothParam <- matrix(NA, ncol = 3, nrow = 30)
for (i in 1:30) {
  res.sc1.smoothParam[i,] <- colMeans(t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[1]]))[,1:3])
}

res.sc2.smoothParam <- matrix(NA, ncol = 3, nrow = 30)
for (i in 1:30) {
  res.sc2.smoothParam[i,] <- colMeans(t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[2]]))[,1:3])
}

res.sc3.smoothParam <- matrix(NA, ncol = 3, nrow = 30)
for (i in 1:30) {
  res.sc3.smoothParam[i,] <- colMeans(t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[3]]))[,1:3])
}

res.sc4.smoothParam <- matrix(NA, ncol = 3, nrow = 30)
for (i in 1:30) {
  res.sc4.smoothParam[i,] <- colMeans(t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[4]]))[,1:3])
}

par(mfrow = c(2,2))
col <- brewer.pal(5, "Dark2")
plot(sqq, rep(NA, 30), ylim = c(0,1), ylab = "Smoothing Parameters", xlab = "Lambda",
     main = "Smoothing parameters with regularisation \n LASSO - Initial values optimised")
abline(h = 0, lwd = 2, col = "grey")
abline(h = 0.4, lwd = 2, col = "red", lty = 6)
lines(sqq,res.sc1.smoothParam[,1], col = col[1], lwd = 2, lty = 2)
lines(sqq,res.sc1.smoothParam[,2], col = col[2], lwd = 2)
lines(sqq,res.sc1.smoothParam[,3], col = col[3], lwd = 2, lty = 4)
legend("topright",
       c("Level smoothing parameter", "Trend smoothing parameter", "Dampening parameter", "DGP Level parameter 0.4"),
       col = c(col[1:3], "red"), lty = c(2, 1, 4, 6), cex = 0.7)

plot(sqq, rep(NA, 30), ylim = c(0,1), ylab = "Smoothing Parameters", xlab = "Lambda",
     main = "Smoothing parameters with regularisation \n LASSO - Initial values backcasting")
abline(h = 0, lwd = 2, col = "grey")
abline(h = 0.4, lwd = 2, col = "red", lty = 6)
lines(sqq,res.sc2.smoothParam[,1], col = col[1], lwd = 2, lty = 2)
lines(sqq,res.sc2.smoothParam[,2], col = col[2], lwd = 2)
lines(sqq,res.sc2.smoothParam[,3], col = col[3], lwd = 2, lty = 4)
legend("topright",
       c("Level smoothing parameter", "Trend smoothing parameter", "Dampening parameter", "DGP Level parameter 0.4"),
       col = c(col[1:3], "red"), lty = c(2, 1, 4, 6), cex = 0.7)

plot(sqq, rep(NA, 30), ylim = c(0,1), ylab = "Smoothing Parameters", xlab = "Lambda",
     main = "Smoothing parameters with regularisation \n Adaptive - Initial values optimised")
abline(h = 0, lwd = 2, col = "grey")
abline(h = 0.4, lwd = 2, col = "red", lty = 6)
lines(sqq,res.sc3.smoothParam[,1], col = col[1], lwd = 2, lty = 2)
lines(sqq,res.sc3.smoothParam[,2], col = col[2], lwd = 2)
lines(sqq,res.sc3.smoothParam[,3], col = col[3], lwd = 2, lty = 4)
legend("topright",
       c("Level smoothing parameter", "Trend smoothing parameter", "Dampening parameter", "DGP Level parameter 0.4"),
       col = c(col[1:3], "red"), lty = c(2, 1, 4, 6), cex = 0.7)

plot(sqq, rep(NA, 30), ylim = c(0,1), ylab = "Smoothing Parameters", xlab = "Lambda",
     main = "Smoothing parameters with regularisation \n Adaptive - Initial values backcasting")
abline(h = 0, lwd = 2, col = "grey")
abline(h = 0.4, lwd = 2, col = "red", lty = 6)
lines(sqq,res.sc4.smoothParam[,1], col = col[1], lwd = 2, lty = 2)
lines(sqq,res.sc4.smoothParam[,2], col = col[2], lwd = 2)
lines(sqq,res.sc4.smoothParam[,3], col = col[3], lwd = 2, lty = 4)
legend("topright",
       c("Level smoothing parameter", "Trend smoothing parameter", "Dampening parameter", "DGP Level parameter 0.4"),
       col = c(col[1:3], "red"), lty = c(2, 1, 4, 6), cex = 0.7)

par(mfrow = c(1,1))


res.sc1.iniValues <- matrix(NA, ncol = 2, nrow = 30)
for (i in 1:30) {
  res.sc1.iniValues[i,] <- colMeans(t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[1]]))[,4:5])
}

res.sc2.iniValues <- matrix(NA, ncol = 2, nrow = 30)
for (i in 1:30) {
  res.sc2.iniValues[i,] <- colMeans(t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[2]]))[,4:5])
}

res.sc3.iniValues <- matrix(NA, ncol = 2, nrow = 30)
for (i in 1:30) {
  res.sc3.iniValues[i,] <- colMeans(t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[3]]))[,4:5])
}

res.sc4.iniValues <- matrix(NA, ncol = 2, nrow = 30)
for (i in 1:30) {
  res.sc4.iniValues[i,] <- colMeans(t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[4]]))[,4:5])
}

plot(sqq, rep(NA, 30), ylim = c(-3,3), ylab = "Initial Values", xlab = "Lambda",
     main = "Initial values with regularisation - LASSO")
abline(h = 0, lwd = 2, col = "grey")
lines(sqq,res.sc1.iniValues[,1], col = col[1], lwd = 2, lty = 3)
lines(sqq,res.sc1.iniValues[,2], col = col[2], lwd = 2, lty = 5)
legend("topright",
       c("Level initial values", "Trend initial values"),
       col = c(col[4:5]), lty = c(3,5), cex = 0.7)

plot(sqq, rep(NA, 30), ylim = c(0,120), ylab = "Initial Values", xlab = "Lambda",
     main = "Initial values with backcasting - LASSO")
abline(h = 0, lwd = 2, col = "grey")
lines(sqq,res.sc2.iniValues[,1], col = col[1], lwd = 2, lty = 3)
lines(sqq,res.sc2.iniValues[,2], col = col[2], lwd = 2, lty = 5)
legend("topright",
       c("Level initial values", "Trend initial values"),
       col = c(col[4:5]), lty = c(3,5), cex = 0.7)

plot(sqq, rep(NA, 30), ylim = c(-3,3), ylab = "Initial Values", xlab = "Lambda",
     main = "Initial values with regularisation - Adaptive")
abline(h = 0, lwd = 2, col = "grey")
lines(sqq,res.sc3.iniValues[,1], col = col[1], lwd = 2, lty = 3)
lines(sqq,res.sc3.iniValues[,2], col = col[2], lwd = 2, lty = 5)
legend("topright",
       c("Level initial values", "Trend initial values"),
       col = c(col[4:5]), lty = c(3,5), cex = 0.7)

plot(sqq, rep(NA, 30), ylim = c(0,120), ylab = "Initial Values", xlab = "Lambda",
     main = "Initial values with backcasting - Adaptive")
abline(h = 0, lwd = 2, col = "grey")
lines(sqq,res.sc4.iniValues[,1], col = col[1], lwd = 2, lty = 3)
lines(sqq,res.sc4.iniValues[,2], col = col[2], lwd = 2, lty = 5)
legend("topright",
       c("Level initial values", "Trend initial values"),
       col = c(col[4:5]), lty = c(3,5), cex = 0.7)


res.sc1.truemodel <- NULL
for (i in 1:30) {
  x <- t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[1]]))[,1:3]
  res.sc1.truemodel <- rbind(res.sc1.truemodel, colSums(apply(x, 2, function(y) ifelse(y>1e-3, 1, 0)))/500)
}

res.sc2.truemodel <- NULL
for (i in 1:30) {
  x <- t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[2]]))[,1:3]
  res.sc2.truemodel <- rbind(res.sc2.truemodel, colSums(apply(x, 2, function(y) ifelse(y>1e-3, 1, 0)))/500)
}

res.sc3.truemodel <- NULL
for (i in 1:30) {
  x <- t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[3]]))[,1:3]
  res.sc3.truemodel <- rbind(res.sc3.truemodel, colSums(apply(x, 2, function(y) ifelse(y>1e-3, 1, 0)))/500)
}

res.sc4.truemodel <- NULL
for (i in 1:30) {
  x <- t(apply(t(sapply(res.exp1, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[4]]))[,1:3]
  res.sc4.truemodel <- rbind(res.sc4.truemodel, colSums(apply(x, 2, function(y) ifelse(y>1e-3, 1, 0)))/500)
}

plot(sqq, res.sc1.truemodel[,1], type = "l", col = col[1], lty = 2, lwd = 2,
     ylab = "Non-zero parameters in percentage", xlab = "Lambda", main = "LASSO with optimal initial values")
lines(sqq, res.sc1.truemodel[,2], type = "l", col = col[2], lty = 1, lwd = 2)
lines(sqq, res.sc1.truemodel[,3], type = "l", col = col[3], lty = 4, lwd = 2)
legend("topright",
       c("Level smoothing parameter", "Trend smoothing parameter", "Dampening parameter"),
       col = c(col[1:3]), lty = c(2, 1, 4), cex = 0.7)

plot(sqq, res.sc2.truemodel[,1], type = "l", col = col[1], lty = 2, lwd = 2,
     ylab = "Non-zero parameters in percentage", xlab = "Lambda", main = "LASSO with backcast initial values")
lines(sqq, res.sc2.truemodel[,2], type = "l", col = col[2], lty = 1, lwd = 2)
lines(sqq, res.sc2.truemodel[,3], type = "l", col = col[3], lty = 4, lwd = 2)
abline(v = 0.3856620421, col = "black")
abline(v = sqq[which(colMeans(res.sc2.outRMSE) == min(colMeans(res.sc2.outRMSE)))], col = "black", lty = 3)
abline(v = sqq[22], col = "black", lty = 5)
legend("topright",
       c("Level smoothing parameter", "Trend smoothing parameter", "Dampening parameter", 
         "Lambda at 0.38", "Optimal Lambda at 0.04", "one-se Lambda at 0.08"),
       col = c(col[1:3], "black"), lty = c(2, 1, 4, 1), cex = 0.7)

plot(sqq, res.sc3.truemodel[,1], type = "l", col = col[1], lty = 2, lwd = 2, ylim = c(0,1),
     ylab = "Non-zero parameters in percentage", xlab = "Lambda", main = "Adaptive with optimal initial values")
lines(sqq, res.sc3.truemodel[,2], type = "l", col = col[2], lty = 1, lwd = 2)
lines(sqq, res.sc3.truemodel[,3], type = "l", col = col[3], lty = 4, lwd = 2)
legend("topright",
       c("Level smoothing parameter", "Trend smoothing parameter", "Dampening parameter"),
       col = c(col[1:3]), lty = c(2, 1, 4), cex = 0.7)

plot(sqq, res.sc4.truemodel[,1], type = "l", col = col[1], lty = 2, lwd = 2, ylim = c(0,1),
     ylab = "Non-zero parameters in percentage", xlab = "Lambda", main = "Adaptive with backcast initial values")
lines(sqq, res.sc4.truemodel[,2], type = "l", col = col[2], lty = 1, lwd = 2)
lines(sqq, res.sc4.truemodel[,3], type = "l", col = col[3], lty = 4, lwd = 2)
abline(v = 0.2807216204 , col = "black")
abline(v = sqq[which(colMeans(res.sc4.outRMSE) == min(colMeans(res.sc4.outRMSE)))], col = "black", lty = 3)
abline(v = sqq[22], col = "black", lty = 5)
legend("topright",
       c("Level smoothing parameter", "Trend smoothing parameter", "Dampening parameter", 
         "Lambda at 0.28", "Optimal Lambda at 0.03", "one-se Lambda at 0.08"),
       col = c(col[1:3], "black"), lty = c(2, 1, 4, 1), cex = 0.7)

rbind(
c("Scenario", "Optimal Lambda", "Validation set RMSE", "one-se Lambda", "Validation set RMSE"),
c("LASSO with backcast initial values", sqq[20], min(colMeans(res.sc2.outRMSE)), sqq[22], colMeans(res.sc2.outRMSE)[22]),
c("Adaptive with backcast initial values", sqq[19], min(colMeans(res.sc4.outRMSE)), sqq[22], colMeans(res.sc4.outRMSE)[22]))

xtable(
rbind(cbind(sqq[25], min(colMeans(res.sc1.outRMSE)), sqq[27], colMeans(res.sc1.outRMSE)[27]),
      cbind(sqq[20], min(colMeans(res.sc2.outRMSE)), sqq[22], colMeans(res.sc2.outRMSE)[22]),
      cbind(sqq[26], min(colMeans(res.sc3.outRMSE)), sqq[28], colMeans(res.sc3.outRMSE)[28]),
      cbind(sqq[20], min(colMeans(res.sc4.outRMSE)), sqq[22], colMeans(res.sc4.outRMSE)[22])), digits =  4)
