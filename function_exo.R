library(smooth)

yData <- sim.es(model = "AAdN", obs = 100, nsim = 1000, frequency = 10, 
                bounds = "usual", persistence = c(0.5, 0.4), phi = 0.7,
                initial = c(100, 2))

sq <- seq(0,1,0.1)
l.sq <- length(sq)

iter <- 50

result.OptLASSO <- vector("list", iter)
result.BcsLASSO <- vector("list", iter)
result.OptPrio <- vector("list", iter)
result.BcsPrio <- vector("list", iter)

start.time <- Sys.time()

iter.run <- function() {
  
  y <- sim.es(model = "AAdN", obs = 100, nsim = 1, frequency = 10, 
              bounds = "usual", persistence = c(0.5, 0.4), phi = 0.7,
              initial = c(100, 2))$data
  
  sq <- c(0, rev(exp(seq.int(from = log(1), to = log(1e-04), by = -(log(1)-log(1e-04))/(100-1)))))
  l.sq <- length(sq)
  
  matLambda.OptLASSO <- matrix(NA, nrow = l.sq, ncol = 17)
  matLambda.BcsLASSO <- matrix(NA, nrow = l.sq, ncol = 17)
  matLambda.OptPrio <- matrix(NA, nrow = l.sq, ncol = 17)
  matLambda.BcsPrio <- matrix(NA, nrow = l.sq, ncol = 17)
  
  for (l in 1:l.sq) {
    
    lambda <<- sq[l]
    
    mat.OptLASSO <- matrix(NA, nrow = 7, ncol = 17)
    mat.BcsLASSO <- matrix(NA, nrow = 7, ncol = 17)
    mat.OptPrio <- matrix(NA, nrow = 7, ncol = 17)
    mat.BcsPrio <- matrix(NA, nrow = 7, ncol = 17)
    
    for (j in 2:8) {
      
      yTrain <- window(y, start = c(1,1), end = c(j,10))
      yHoldout <- window(y, start = c(j+1,1), end = c(j+1,10))
      
      fit.OptLASSO <- adam(yTrain, model = "AAdA", initial = "optimal", loss = "LASSO", 
                           lags = 10, lambda = lambda, h = 10)
      fit.BcsLASSO <- adam(yTrain, model = "AAdA", initial = "backcasting", loss = "LASSO", 
                           lags = 10, lambda = lambda, h = 10)
      fit.OptPrio <- adam(yTrain, model = "AAdA", initial = "optimal", loss = myLoss_OptimalInitials, 
                          lags = 10, h = 10)
      fit.BcsPrio <- adam(yTrain, model = "AAdA", initial = "backcasting", loss = myLoss_BackcastingInitials, 
                          lags = 10, h = 10)
      
      holdOutRMSE.OptLASSO <- sqrt(mean((yHoldout-fit.OptLASSO$forecast)^2))
      holdOutRMSE.BcsLASSO <- sqrt(mean((yHoldout-fit.BcsLASSO$forecast)^2))
      holdOutRMSE.OptPrio <- sqrt(mean((yHoldout-fit.OptPrio$forecast)^2))
      holdOutRMSE.BcsPrio <- sqrt(mean((yHoldout-fit.BcsPrio$forecast)^2))
      
      param.OptLASSO <- c(fit.OptLASSO$persistence, fit.OptLASSO$phi, unlist(fit.OptLASSO$initial))
      param.BcsLASSO <- c(fit.BcsLASSO$persistence, fit.BcsLASSO$phi, unlist(fit.BcsLASSO$initial))
      param.OptPrio <- c(fit.OptPrio$persistence, fit.OptPrio$phi, unlist(fit.OptPrio$initial))
      param.BcsPrio <- c(fit.BcsPrio$persistence, fit.BcsPrio$phi, unlist(fit.BcsPrio$initial))
      
      mat.OptLASSO[j-1,] <- c(holdOutRMSE.OptLASSO, param.OptLASSO)
      mat.BcsLASSO[j-1,] <- c(holdOutRMSE.BcsLASSO, param.BcsLASSO)
      mat.OptPrio[j-1,] <- c(holdOutRMSE.OptPrio, param.OptPrio)
      mat.BcsPrio[j-1,] <- c(holdOutRMSE.BcsPrio, param.BcsPrio)
      
    }
    
    matLambda.OptLASSO[l,] <- colMeans(mat.OptLASSO)
    matLambda.BcsLASSO[l,] <- colMeans(mat.BcsLASSO)
    matLambda.OptPrio[l,] <- colMeans(mat.OptPrio)
    matLambda.BcsPrio[l,] <- colMeans(mat.BcsPrio)

  }
  
  return(list(matLambda.OptLASSO,
              matLambda.BcsLASSO,
              matLambda.OptPrio,
              matLambda.BcsPrio))
  
}

y <- yData$data[,1]
 
system.time({xyz <- iter.run()})
