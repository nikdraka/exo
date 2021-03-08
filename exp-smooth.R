library(smooth)

y <- sim.es(model = "ANA", obs = 100, nsim = 1, frequency = 10, bounds = "usual", persistence = 0.4)
ts.plot(y$data)

mat.fit <- array(NA, c(length(seq(0,1,0.1)), 3))
for (i in seq(0, 1, 0.1)) {
    fit <- adam(y, model = "AAdA", initial = "optimal", loss = "LASSO", lambda = i)
    mat.fit[i*10+1,] <- c(fit$persistence, fit$phi)
}


myLoss_OptimalInitials <- function(actual, fitted, B) {
    
    # Make sure that we use seasonal models
    
    nB <- length(B)
    obsActual <- length(actual)
    yDenominator <- max(sd(diff(actual)),1)
    
    priority <- c(c(2, 1, 0, 1, 2, 1), c(rep(0, nB-6)))
    B <- B * priority

    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)

    CFValue <- (1-lambda)*loss + lambda*sum(abs(B))

    return(CFValue)
}

myLoss_BackcastingInitials <- function(actual, fitted, B) {

    obsActual <- length(actual)
    yDenominator <- max(sd(diff(actual)),1)

    priority <- c(2, 1, 0, 1)
    B <- B * priority

    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)

    CFValue <- (1-lambda)*loss + lambda*sum(abs(B))

    return(CFValue)
}

for (lambda in seq(0,1,0.1)) {

    lambda <<- lambda
    print(adam(y, model = "AAdA", initial = "optimal", loss = myLoss_OptimalInitials, lags = 10)$lossValue)
    print(adam(y, model = "AAdA", initial = "backcasting", loss = myLoss_BackcastingInitials, lags = 10)$lossValue)

}

lambda <- 0.1
fit <- adam(y, model = "AAdA", initial = "optimal", loss = "LASSO", lags = 10, lambda = lambda)
fit <- adam(y, model = "AAdA", initial = "backcasting", loss = "LASSO", lags = 7, lambda = lambda)
fit <- adam(y, model = "AAdA", initial = "optimal", loss = myLoss_OptimalInitials, lags = 10)
fit <- adam(y, model = "AAdA", initial = "backcasting", loss = myLoss_BackcastingInitials, lags = 10, h = 10)
