nObs <- 70
y.exp <- sim.es(model = "ANN", obs = nObs, nsim = 500, frequency = 10, 
                persistence = 0.4, initial = 100,
                bounds = "usual")$data

y.exp.out <- array(NA, c(10,500))
for (i in 1:500) {
  y.exp.out[,i] <- window(y.exp[,i], start = c(7,1))
}

y.exp.ins <- ts(matrix(NA, ncol = 500, nrow = 60), frequency = 10, start = c(1,1))
for (i in 1:500) {
  y.exp.ins[,i] <- ts(window(y.exp[,i], end = c(6,10)), frequency = 10, start = c(1,1))
}
