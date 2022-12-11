library(mlbench)
library(vcmasf)
library(splines2)
library(pracma)
data("BostonHousing2")

## Analyze Boston Housing data
y = BostonHousing2$medv
y = log(y)
u = BostonHousing2$lstat
u = rank(u) / (length(u) + 1)
X = BostonHousing2[, c('crim', 'rm', 'ptratio', 'nox', 'tax', 'age')]
for (i in 1:6) {
  X[,i] = qnorm(rank(X[,i]) / (length(u) + 1))
  fit = lm(X[,i]~u)
  X[,i] = fit$residuals
  X[,i] = X[,i] / sqrt(mean(X[,i]^2))
}
X = cbind(1, X)
features = c("int", 'crim', 'rm', 'ptratio', 'nox', 'tax', 'age')
degree = 3
fit1 = vcm.asf.predictor.specific(X, y, u, degree=degree, boundary=c(0, 1))

## Get confidence interval for the varying coefficients
n = length(y)
bspline = matrix(0, nrow=n, ncol=0)
for (j in 1:7) {
  uspline = bSpline(u, knots=fit1$knots[[j]], degree=degree, Boundary.knots=c(0, 1), intercept=TRUE)
  bspline = cbind(bspline, X[,j] * uspline)
}
lm_fit = lm(y~bspline-1)
coef = lm_fit$coefficient
sigma = sqrt(mean(lm_fit$residuals^2))
Sigma = inv(t(bspline) %*% bspline)

start = 1
end = 0
sds = matrix(0, ncol=7, nrow=n)
for (i in 1:7) {
  end = start + degree + length(fit1$knots[[i]])
  uspline = bSpline(u, knots=fit1$knots[[i]], degree=degree, Boundary.knots=c(0, 1), intercept=TRUE)
  sds[, i] = sigma * sqrt(diag(uspline %*% Sigma[c(start:end), c(start:end)] %*% t(uspline)))
  start = end + 1
}

# Visualize the confidence interval
par(mfrow=c(2,4), mar=c(4,2,2,1.5))
for (i in 1:7) {
  ind = order(u)
  coef = fit1$coef(u)[ind, i]
  plot(u[ind], coef, main=features[i], type='l', xlab='transformed lstat', ylab='beta',
       ylim=c(min(coef - 2 * sds[ind, i]), max(coef + 2 * sds[ind, i])), lwd=2)
  lines(u[ind], coef - 2 * sds[ind, i], lty=3, lwd=2)
  lines(u[ind], coef + 2 * sds[ind, i], lty=3, lwd=2)
  lines(u[ind], rep(0, n), col='black', lty=2, lwd=2)
}

## 10-fold cross-validation
n = length(u)
nfolds = 10
folds = c(1:n) %% nfolds
f0 = rep(0, n)
f1 = rep(0, n)
for (i in 0:9) {
  training = which(folds != i)
  testing = which(folds == i)
  fit0 = lm(y[training] ~ X[training, 2] + X[training, 3] + X[training, 4] + X[training, 5] + X[training, 6] + X[training, 7] + u[training])
  fit1 = vcm.asf.predictor.specific(X[training, ], y[training], u[training], degree=degree, boundary=c(0, 1))
  for (k in 1:7) {
    f0[testing] = f0[testing] + X[testing, k] * fit0$coefficients[k]
  }
  f0[testing] = f0[testing] + u[testing] * fit0$coefficients[8]
  f1[testing] = fit1$predict(X[testing, ], u[testing])
}

## Compare MSE between simple linear model and predictor-specific varying coefficient model
print(c(mean((exp(y)-exp(f0))^2), mean((exp(y)-exp(f1))^2)))
