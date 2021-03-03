
createInit <- function(y, type = c("optimized", "backcasting")) {
  
  freq.ts <- frequency(y)
  
  alpha <- 0.1
  beta <- 0.05
  phi <- 0.95
  
  if (type == "optimized") {
    x0 <- unname(coef(lm(head(y,freq.ts)~seq(1,freq.ts,1))))
    x0[1] <- (x0[1] - mean(head(y, freq.ts)))/sd(head(y, freq.ts))
    x0[2] <- x0[2]/sd(head(y, freq.ts))
    ub.x0 <- rep(3, 2)
    lb.x0 <- rep(-3, 2)
  } else if (type == "backcasting") {
    x0 <- NULL
    ub.x0 <- NULL
    lb.x0 <- NULL
  }
  
  return(list(initials = c(alpha, beta, phi, x0),
              lb = c(0,0,0,lb.x0),
              ub = c(1,1,1,ub.x0)))
  
}

x0 <- createInit(y, type = "optimized")

lossFun <- function(x0, y, type = c("optimized", "backcasting"), 
                    penalty = c("lasso", "adaptive", "enet", "adaptive-enet"), 
                    lambda = 0, elastic = NULL) {
  
  nTrain <- length(y)
  
  if (x0[1] > 1 || x0[2] > 1) {
    loss <- 1e+300
  } else if (x0[2] > x0[1]) {
    loss <- 1e+300
  }

  if (x0[3] & (x0[2] < 0 || x0[2] > 1)) {
    loss <- 1e+300
  }
  
  matH <- matrix(c(1,1), ncol = 2, byrow = TRUE)
  matF <- matrix(c(1, x0[3],
                   0, x0[3]), ncol = 2, byrow = TRUE)
  matG <- matrix(c(x0[1], x0[2]), ncol = 1, byrow = TRUE)
  
  fitted <- matrix(NA, ncol = 1, nrow = nTrain)
  statevec <- matrix(NA, ncol = 2, nrow = nTrain)
  
  if (type == "optimized") {
    
    s0 <- matrix(c(x0[4], x0[5]), ncol = 1, byrow = TRUE)
    
  } else if (type == "backcasting") {
    
    freq.ts <- frequency(y)
    s0 <- matrix(unname(coef(lm(head(y,freq.ts)~seq(1,freq.ts,1)))), ncol = 1, byrow = TRUE)
    
  }
  
  if (type == "optimized") {
    l0 <- s0
    for (i in 1:nTrain) {
      
      yhat <- matH %*% (l0)
      ehat <- y[i] - yhat
      
      l1 <- matF %*% (s0) + matG %*% ehat
      
      fitted[i] <- yhat
      statevec[i,] <- t(l1)
      l0 <- l1
    }
    
  } else if (type == "backcasting") {
    l0 <- s0
    for (i in 1:nTrain) {
      
      yhat <- matH %*% l0
      ehat <- y[i] - yhat
      
      l1 <- matF %*% l0 + matG %*% ehat
      
      fitted[i] <- yhat
      statevec[i,] <- l1
      l0 <- rbind(statevec[i,1], x0[3]*statevec[i,2])
    }
    
  }
  
  error <- y - fitted
  mse <- (1-lambda)*sqrt(mean(error^2))/var(diff(y))
  
  if (type == "optimized") {
    est.param <- c(matG, matF[1,2], s0)
    w <- 0.5*c(2, 1, 0, 3, 3)
  } else if (type == "backcasting") {
    est.param <- c(matG, matF[1,2])
    w <- 0.5^c(2,1,0)
  }
  
  if (penalty == "lasso") {
    
    shrink.penalty <- lambda*sum(abs(est.param))
    
  } else if (penalty == "adaptive") {
    
    # weight.model <- es(y, model = "AAdN", loss = "MSE")
    # if (type == "optimized") {
    #   weight.model$initial[1] <- (weight.model$initial[1] - mean(head(y, freq.ts)))/sd(head(y, freq.ts))
    #   weight.model$initial[2] <- (weight.model$initial[2])/sd(head(y, freq.ts))
    #   w <- 1/abs(c(weight.model$persistence, weight.model$phi, weight.model$initial))
    # } else if (type == "backcasting") {
    #   est.param <- c(matG, matF[1,2])
    #   weight.model <- es(y, model = "AAdN", loss = "MSE")
    #   w <- 1/abs(c(weight.model$persistence, weight.model$phi))
    # }
    # 
    # w.estparam <- rbind(is.finite(w), w, est.param)
    # w.estparam <- w.estparam[,-which(apply(w.estparam, 1, function(y) y== FALSE))]
    
    shrink.penalty <- lambda*(w %*% est.param)
    
  } else if (penalty == "enet") {
    
    shrink.penalty <- elastic*lambda*sum(abs(est.param)) + (1-elastic)*lambda*sum((est.param)^2)
    
  } else if (penalty == "adaptive-enet") {
    
    weight.model <- es(y, model = "AAdN", loss = "MSE")
    if (type == "optimized") {
      weight.model$initial[1] <- (weight.model$initial[1] - mean(head(y, freq.ts)))/sd(head(y, freq.ts))
      weight.model$initial[2] <- (weight.model$initial[2])/sd(head(y, freq.ts))
      w <- 1/abs(c(weight.model$persistence, weight.model$phi, weight.model$initial))
    } else if (type == "backcasting") {
      w <- 1/abs(c(weight.model$persistence, weight.model$phi))
    }
    
    w.estparam <- rbind(is.finite(w), w, est.param)
    w.estparam <- w.estparam[2:3,-which(apply(w.estparam, 1, function(y) y== FALSE))]
    shrink.penalty <- elastic*lambda*sum(abs(sum(apply(w.estparam, 2, prod)))) + (1-elastic)*lambda*sum((sum(apply(w.estparam, 2, prod)))^2)
    
  }
  
  loss <- mse + shrink.penalty
  return(c(loss))
  
}

lossFun(x0$initials, y, type = "optimized", penalty = "adaptive", lambda = 0.5, elastic = NULL)
lossFun(x0$initials, y, type = "optimized", penalty = "lasso", lambda = 0.5, elastic = NULL)

y <- sim.es(model = "ANN", obs = 70, nsim = 1, frequency = 10, 
            persistence = 0.4, initial = 100,
            bounds = "usual")$data



loopLambdaCV <- function(data = y, type = c("optimized", "backcasting"), 
                         penalty = c("lasso", "adaptive", "enet", "adaptive-enet"),
                         sq = seq(0, 1, 0.1)) {
  
  y <- data
  type <- type
  penalty <- penalty
  
  sq <- sq
  # mat.lambda <- matrix(NA, ncol = 9, nrow = length(sq))
  mat.lambda <- NULL
  for (lamb in 1:length(sq)) {
    
    # mat.origin <- matrix(NA, ncol = 8, nrow = 6)
    mat.origin <- NULL
    
    for (origin in 1:6) {
      yTrain <- window(y, end = c(origin,10))
      yVal <- window(y, start = c(origin+1,1), end = c(origin+1,10))
      nTrain <- length(yTrain)
      
      opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel" = 1e-08, "maxeval" = 1000)
      
      x0 <- createInit(yTrain, type = type)
      
      fit <- nloptr(x0$initials, function(x) lossFun(x, yTrain, type = type, 
                                                     penalty = penalty, lambda = sq[lamb], elastic = NULL),
                    lb = x0$lb, ub = x0$ub, opts = opts)
      
      mat <- matrixMaker(fit$solution, yTrain, type = type)
      
      freq.ts <- frequency(yTrain)
      fitted <- matrix(NA, ncol = 1, nrow = nTrain+1)
      statevec <- matrix(NA, ncol = 2, nrow = nTrain+1)
      
      if (type == "optimized") {
        mat$initialValue[1] <- mat$initialValue[1] * sd(head(yTrain, freq.ts)) + mean(head(yTrain, freq.ts))
        mat$initialValue[2] <- mat$initialValue[2] * sd(head(yTrain, freq.ts))
      }
      
      if (type == "optimized") {
        l0 <- mat$initialValue
        for (i in 1:(nTrain+1)) {
          
          yhat <-mat$measurementMatrix %*% l0
          ehat <- yTrain[i] - yhat
          
          l1 <- mat$transitionMatrix %*% l0 + mat$persistenceMatrix %*% ehat
          
          fitted[i] <- yhat
          statevec[i,] <- t(l1)
          l0 <- l1
        }
        
      } else if (type == "backcasting") {
        l0 <- mat$initialValue
        for (i in 1:(nTrain+1)) {
          
          yhat <- mat$measurementMatrix %*% l0
          ehat <- y[i] - yhat
          
          l1 <- mat$transitionMatrix %*% l0 + mat$persistenceMatrix %*% ehat
          
          fitted[i] <- yhat
          statevec[i,] <- l1
          l0 <- rbind(statevec[i,1], fit$solution[3]*statevec[i,2])
        }
      }
      
      fcst <- gen.forecast(fit$solution, statevec, h = 10)
      fcst.error <- yVal - fcst
      outRMSE <- sqrt(mean((fcst.error)^2))
      
      mat.origin <- rbind(mat.origin, c(fit$solution, mat$initialValue, outRMSE))
      # mat.origin[origin,] <- c(fit$solution, mat$initialValue, outRMSE)
    }
    
    # print(c(apply(mat.origin, 2, mean)))
    
    mat.lambda <- rbind(mat.lambda, c(apply(mat.origin, 2, mean),
                                      sd(mat.origin[,ncol(mat.origin)])/sqrt(6)))
    # mat.lambda[lamb,] <- c(apply(mat.origin, 2, mean), sd(mat.origin[,6])/sqrt(6))
  }
  return(mat.lambda)
}

fun.iter <- function(sq = seq(0, 1, 0.1)) {
  
  sq <- sq
  
  y <- sim.es(model = "ANN", obs = 80, nsim = 1, frequency = 10, 
              persistence = 0.4, initial = 100,
              bounds = "usual")$data
  
  coba.lasso.optimal <- loopLambdaCV(y, type = "optimized", penalty = "lasso", sq = sq)
  coba.lasso.backcast <- loopLambdaCV(y, type = "backcasting", penalty = "lasso", sq = sq)
  
  coba.adapt.optimal <- loopLambdaCV(y, type = "optimized", penalty = "adaptive", sq = sq)
  coba.adapt.backcast <- loopLambdaCV(y, type = "backcasting", penalty = "adaptive", sq = sq)
  
  return(list(coba.lasso.optimal,
              coba.lasso.backcast,
              coba.adapt.optimal,
              coba.adapt.backcast))
  
}

sqq <- rev(exp(seq.int(from = log(1), to = log(1e-04), by = -(log(1)-log(1e-04))/(10-1))))
system.time({coba.fun.iter <- fun.iter(sq = sqq)})
