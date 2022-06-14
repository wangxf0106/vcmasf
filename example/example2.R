library(vcmasf)
library(mvtnorm)
library(rdist)

example2 = function(n=50, ni=30, skip=0.6) 
{
  y = c()
  f = c()
  u = c()
  X = matrix(0, ncol=500, nrow=0)
  B = matrix(0, ncol=6, nrow=0)
  st = seq(0, ni-1, 1)
  for (i in 1:n) {
    ot = st + runif(ni, 0, 1)
    not_skip = rbinom(ni, 1, 1-skip)
    if (sum(not_skip) == 0)
      next
    ot = ot[not_skip == 1]
    m = length(ot)
    x = matrix(0, ncol=500, nrow=m)
    x[,1] = (ot+0.5)/10 + runif(m, 0, 2)
    for (j in 2:5)
      x[,j] = rnorm(m) * sqrt((1+x[,1])/(2+x[,1]))
    x[,6] = rnorm(m) + 3*exp((ot+0.5)/30)
    cov = 4*exp(-as.matrix(rdist(ot)))
    x[,7:500] = t(rmvnorm(494, rep(0,m), cov))
    z = as.numeric(rmvnorm(1, rep(0,m), cov))
    e = rnorm(m, 0, 2)
    
    b = matrix(0, ncol=6, nrow=m)
    b[,1] = 15+20*sin(pi*(ot+0.5)/15)
    b[,2] = 15+20*cos(pi*(ot+0.5)/15)
    b[,3] = 2-3*sin(pi*(ot-24.5)/15)
    b[,4] = 2-3*cos(pi*(ot-24.5)/15)
    b[,5] = 6-0.2*(ot+0.5)^2
    b[,6] = -4+(19.5-ot)^3/2000
    f = c(f, rowSums(b*x[,1:6]))
    y = c(y, rowSums(b*x[,1:6]) + z + e)
    u = c(u, ot)
    X = rbind(X, x)
    B = rbind(B, b)
  }
  return(list(y=y, X=X, u=u, f=f, B=B))
}

data = example2()
boundary = c(0, 30)
res = vcm.variable.selection(data$X, data$y, data$u, boundary=boundary)
print(res$ind)
print(cor(data$y, res$predict(data$X, data$u)))

ntotal = length(res$ind)
par(mfrow=c(ceiling(ntotal/3), 3))
ind = order(data$u)
B = res$coef(data$u)
for (i in 1:ntotal) {
  plot(data$u[ind], data$B[ind, res$ind[i]], type='l', xlab='u', ylab=paste0('B', i))
  lines(data$u[ind], B[ind, i], col='red')
}