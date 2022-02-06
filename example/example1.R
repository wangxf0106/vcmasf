library(vcmasf)
library(mvtnorm)
library(rdist)

example1 = function(n=200, ni=20, skip=0.6) {
  f = c()
  y = c()
  u = c()
  X = matrix(0, ncol=4, nrow=0)
  B = matrix(0, ncol=4, nrow=0)
  st = seq(0, ni-1, 1)
  for (i in 1:n) {
    ot = st + runif(ni)
    not_skip = rbinom(ni, 1, 1-skip)
    if (sum(not_skip) == 0)
      next
    ot = ot[not_skip == 1]
    m = length(ot)
    x = matrix(0, ncol=4, nrow=m)
    x[,1] = 1.0
    x[,2] = rbinom(m, 1, 0.6)
    x[,3] = runif(m, 0, 2) + ot/10
    x[,4] = rnorm(m) * sqrt((1+x[,3]/(2+x[,3])))
    cov = 4*exp(-as.matrix(rdist(ot)))
    v = as.numeric(rmvnorm(1, rep(0,m), cov))
    e = rnorm(m, 0, 2)
    b = matrix(0, ncol=4, nrow=m)
    b[,1] = 1+3.5*sin(ot-3)
    b[,2] = 2-5*cos((3*ot-1)/4)
    b[,3] = 4 - ((ot-12)/5)^2
    b[,4] = 1+ot/8+4.6*(1-ot/10)^3
    f = c(f, rowSums(b*x))
    y = c(y, rowSums(b*x) + v + e)
    u = c(u, ot)
    X = rbind(X, x)
    B = rbind(B, b)
  }
  return(list(y=y, X=X, u=u, B=B, f=f))
}
  
data = example1()
boundary = c(0, 20)
fit1 = vcm.asf.onestep(data$X, data$y, data$u, boundary=boundary)
fit2 = vcm.asf.twostep(data$X, data$y, data$u, boundary=boundary)
fit3 = vcm.asf.equidistant(data$X, data$y, data$u, boundary=boundary)
B1 = fit1$coef(data$u)
B2 = fit2$coef(data$u)
B3 = fit3$coef(data$u)

par(mfrow=c(2, 2))
for (i in 1:4) {
  plot(data$u[ind], data$B[ind, i], type='l', xlab='u', ylab=paste0('B', i))
  lines(data$u[ind], B1[ind, i], col='red')
  lines(data$u[ind], B2[ind, i], col='blue')
  lines(data$u[ind], B3[ind, i], col='green')
}
