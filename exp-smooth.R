library(smooth)

y <- sim.es(model = "ANA", obs = 100, nsim = 1, frequency = 10, bounds = "usual", persistence = 0.4)
ts.plot(y$data)

mat.fit <- array(NA, c(length(seq(0,1,0.1)), 3))
for (i in seq(0, 1, 0.1)) {
    fit <- adam(y, model = "AAdA", initial = "optimal", loss = "LASSO", lambda = i)
    mat.fit[i*10+1,] <- c(fit$persistence, fit$phi)
}


myLoss_OptimalInitials <- function(actual, fitted, B) {

    obsActual <- length(actual)
    yDenominator <- max(sd(diff(actual)),1)

    if (length(B) == 5) {
        priority <- 0.5^c(1,0,0,1,0)
        B <- B * priority
    } else if (length(B) == 7) {
        priority <- 0.5^c(2,1,0,1,2,1,0)
        B <- B * priority
    }

    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)

    CFValue <- (1-lambda)*loss + lambda*sum(abs(B))

    return(CFValue)
}

myLoss_BackcastingInitials <- function(actual, fitted, B) {

    obsActual <- length(actual)
    yDenominator <- max(sd(diff(actual)),1)

    if (length(B) == 3) {
        priority <- 0.5^c(1,0,0)
        B <- B * priority
    } else if (length(B) == 5) {
        priority <- 0.5^c(2,1,0,1)
        B <- B * priority
    }

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


adam(y, model = "AAdA", initial = "backcasting", loss = "LASSO", lags = 10, lambda = lambda)
adam(y, model = "AAdA", initial = "optimal", loss = myLoss_OptimalInitials, lags = 10)
fit <- adam(y, model = "AAdA", initial = "backcasting", loss = myLoss_BackcastingInitials, lags = 10)
