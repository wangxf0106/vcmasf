library(vcmasf)
library(mvtnorm)
library(rdist)
library(tvReg)

## Generate simulation data for Section 4.1
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

## Visualize the fitted coefficients by equidistant and predictor-specific method
data = example1()
boundary = c(0, 20)
fit1 = vcm.asf.predictor.specific(data$X, data$y, data$u, boundary=boundary)
fit2 = vcm.asf.equidistant(data$X, data$y, data$u, boundary=boundary)
B1 = fit1$coef(data$u)
B2 = fit2$coef(data$u)
png('example1_coefficient_20221210.png', width=460, height=115)
par(mfrow=c(1,4), mar=c(4,2,2,1.5), ps=16, font.main=16)
for (i in 1:4) {
  ind = order(data$u)
  lower = min(B1[,i], B2[,i]) - 1
  upper = max(B1[,i], B2[,i]) + 1
  plot(data$u[ind], data$B[ind, i], type='l', xlab='u', ylab=paste0('B', i), ylim=c(lower, upper), main=paste0('beta', i))
  grid(6, NA, col='grey')
  lines(data$u[ind], B1[ind, i], lty=3)
  lines(data$u[ind], B2[ind, i], lty=2)
  points(fit1$knots[[i]], rep(upper - 0.5, length(fit1$knots[[i]])), pch=8)
  points(fit2$knots, rep(lower + 0.5, length(fit2$knots)), pch=6)
}
dev.off()

## Compare MSE for coefficients with 1000 repititions
nrep = 1000
bmse1 = matrix(0, nrow=nrep, ncol=4)
bmse2 = matrix(0, nrow=nrep, ncol=4)
bmse3 = matrix(0, nrow=nrep, ncol=4)
bmse4 = matrix(0, nrow=nrep, ncol=4)
nknot1 = rep(0, nrep)
nknot2 = rep(0, nrep)
nknot3 = matrix(0, nrow=nrep, ncol=4)

for (i in 1:nrep)
{
  data = example1()
  boundary = c(0, 20)
  data_tvlm = data.frame(x1 = data$X[,1], x2 = data$X[,2], x3=data$X[,3], x4=data$X[,4], y=data$y)
  coef_tvlm = coef(res_tvlm)
  fit1 = vcm.asf.equidistant(data$X, data$y, data$u, boundary=boundary)
  fit2 = vcm.asf.global(data$X, data$y, data$u, boundary=boundary)
  fit3 = vcm.asf.predictor.specific(data$X, data$y, data$u, boundary=boundary)
  fit4 = tvLM(y ~ 0 + x1 + x2 + x3 + x4, z=data$u, data=data_tvlm)
  B1 = fit1$coef(data$u)
  B2 = fit2$coef(data$u)
  B3 = fit3$coef(data$u)
  bmse1[i, ] = colMeans((data$B - fit1$coef(data$u))^2)
  bmse2[i, ] = colMeans((data$B - fit2$coef(data$u))^2)
  bmse3[i, ] = colMeans((data$B - fit3$coef(data$u))^2)
  bmse4[i, ] = colMeans((data$B - fit4$coefficients)^2)
  nknot1[i] = length(fit1$knots)
  nknot2[i] = length(fit2$knots)
  nknot3[i, 1] = length(fit3$knots[[1]])
  nknot3[i, 2] = length(fit3$knots[[2]])
  nknot3[i, 3] = length(fit3$knots[[3]])
  nknot3[i, 4] = length(fit3$knots[[4]])
}

brg = c(7, 10, 5.76, 6.7)
## MSE for coefficient
## Equidistant spline fitting
print(colMeans(bmse1) / brg^2 * 100)
print(sqrt(colMeans(bmse1^2) - colMeans(bmse1)^2) / brg^2 * 100)
## Global adaptive spline fitting
print(colMeans(bmse2) / brg^2 * 100)
print(sqrt(colMeans(bmse2^2) - colMeans(bmse2)^2) / brg^2 * 100)
## Predictor-specific adaptive spline fitting
print(colMeans(bmse3) / brg^2 * 100)
print(sqrt(colMeans(bmse3^2) - colMeans(bmse3)^2) / brg^2 * 100)
## Kernel fitting
print(colMeans(bmse4) / brg^2 * 100)
print(sqrt(colMeans(bmse4^2) - colMeans(bmse4)^2) / brg^2 * 100)
## Number of knots selected
## Equidistant spline fitting
print(mean(nknot1))
## Global adaptive spline fitting
print(mean(nknot2))
## Predictor-specific adaptive spline fitting
print(colMeans(nknot3))
