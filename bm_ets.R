benchmark.ets <- function(y) {
  mat.acc <- matrix(NA, ncol = 1, nrow = 6)
  for (origin in 1:6) {
    yTrain <- window(y, end = c(origin,10))
    yVal <- window(y, start = c(origin+1,1), end = c(origin+1,10))
    nTrain <- length(yTrain)
    
    # f <- forecast(ets(yTrain, additive.only = TRUE), h = 10)$mean
    f <- es(yTrain, model = "XXX", h = 10)$forecast
    ef <- yVal - f
    acc.rmse <- sqrt(mean(ef^2))
    mat.acc[origin,] <- acc.rmse
  }
  return(mean(mat.acc))
}

system.time({benchmark.ets(y)})
