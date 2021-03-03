res.sc1.outRMSE.bm <- matrix(NA, ncol = 30, nrow = 500)
for (i in 1:30) {
  x <- t(apply(t(sapply(res.exp1.bm, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[1]]))[,9]
  y <- t(apply(t(sapply(res.exp1.bm, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[1]]))[,1]
  res.sc1.outRMSE.bm[,i] <- x/y
}

res.sc2.outRMSE.bm <- matrix(NA, ncol = 30, nrow = 500)
for (i in 1:30) {
  x <- t(apply(t(sapply(res.exp1.bm, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[2]]))[,8]
  y <- t(apply(t(sapply(res.exp1.bm, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[2]]))[,1]
  res.sc2.outRMSE.bm[,i] <- x/y
}

res.sc3.outRMSE.bm <- matrix(NA, ncol = 30, nrow = 500)
for (i in 1:30) {
  x <- t(apply(t(sapply(res.exp1.bm, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[3]]))[,9]
  y <- t(apply(t(sapply(res.exp1.bm, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[3]]))[,1]
  res.sc3.outRMSE.bm[,i] <- x/y
}

res.sc4.outRMSE.bm <- matrix(NA, ncol = 30, nrow = 500)
for (i in 1:30) {
  x <- t(apply(t(sapply(res.exp1.bm, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[4]]))[,8]
  y <- t(apply(t(sapply(res.exp1.bm, function(x) sapply(x, function(y) y[i,]))), 1, function(z) z[[4]]))[,1]
  res.sc4.outRMSE.bm[,i] <- x/y
}

par(mfcol=c(2,2))
boxplot(res.sc1.outRMSE.bm, main = "Rel-OutRMSE in Validation Set \n LASSO with optimal initial values",
        ylab = "Out-RMSE", xlab = "Lambda", xaxt = "n", ylim = c(0,20))
lines(colMeans(res.sc1.outRMSE.bm), col = "red", type = "o", pch = 19, bg = "red")
axis(1,at=1:30,labels=round(sqq, 4),las=2, cex.axis = 0.8)
abline(h = 1, col = "grey")

boxplot(res.sc2.outRMSE.bm, main = "Rel-OutRMSE in Validation Set \n LASSO with backcast initial values",
        ylab = "Out-RMSE", xlab = "Lambda", xaxt = "n", ylim = c(0,1.5))
lines(colMeans(res.sc2.outRMSE.bm), col = "red", type = "o", pch = 19, bg = "red")
axis(1,at=1:30,labels=round(sqq, 4),las=2, cex.axis = 0.8)
abline(h = 1, col = "grey")

boxplot(res.sc3.outRMSE.bm, main = "Rel-OutRMSE in Validation Set \n Adaptive with optimal initial values",
        ylab = "Out-RMSE", xlab = "Lambda", xaxt = "n", ylim = c(0, 100))
lines(colMeans(res.sc3.outRMSE.bm), col = "red", type = "o", pch = 19, bg = "red")
abline(h = 1, col = "grey")
axis(1,at=1:30,labels=round(sqq, 4),las=2, cex.axis = 0.8)

boxplot(res.sc4.outRMSE.bm, main = "Rel-OutRMSE - Adaptive with backcast initial values",
        ylab = "Out-RMSE", xlab = "Lambda", xaxt = "n", ylim = c(0,1.5))
lines(colMeans(res.sc4.outRMSE.bm), col = "red", type = "o", pch = 19, bg = "red")
abline(h = 1, col = "grey")
axis(1,at=1:30,labels=round(sqq,4),las=2, cex.axis = 0.8)

par(mfcol=c(1,1))

boxplot(res.sc2.outRMSE.bm, main = "Rel-OutRMSE in Validation Set \n LASSO with backcast initial values",
        ylab = "Out-RMSE", xlab = "Lambda", xaxt = "n", ylim = c(0,0.5))
lines(colMeans(res.sc2.outRMSE.bm), col = "red", type = "o", pch = 19, bg = "red")
axis(1,at=1:30,labels=round(sqq, 4),las=2, cex.axis = 0.8)

boxplot(res.sc4.outRMSE.bm, main = "Rel-OutRMSE - Adaptive with backcast initial values",
        ylab = "Out-RMSE", xlab = "Lambda", xaxt = "n", ylim = c(0,0.5))
lines(colMeans(res.sc4.outRMSE.bm), col = "red", type = "o", pch = 19, bg = "red")
abline(h = 1, col = "grey")
axis(1,at=1:30,labels=round(sqq,4),las=2, cex.axis = 0.8)