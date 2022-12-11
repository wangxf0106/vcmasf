library(vcmasf)

## Analyze COVID-19 case numbers and environmental factors data
## Data preparation 
covid.environment[which(covid.environment$case_avg < 1), 'case_avg'] = 1
covid.environment['log_case_avg'] = log(covid.environment['case_avg'])
covid.environment['time'] = as.Date(covid.environment$time, format='%Y-%m-%d')
covid.environment['year'] = unclass(as.POSIXlt(covid.environment$time))$year
covid.environment['year'] = covid.environment['year'] - 120
covid.environment['yday'] = unclass(as.POSIXlt(covid.environment$time))$yday 
covid.environment[which(covid.environment$year == 0), 'yday'] = covid.environment[which(covid.environment$year == 0), 'yday'] / 366
covid.environment[which(covid.environment$year != 0), 'yday'] = covid.environment[which(covid.environment$year != 0), 'yday'] / 365
covid.environment['t'] = covid.environment['year'] + covid.environment['yday']

## Visualization of New York County data
par(mfrow=c(2,4), mar=c(4,2,2,1.5))
ny = covid.environment[covid.environment$county == "New York", ]
cols = c("temp", "dew", "wind", "precipitation", "humidity", "pm25", "ozone")
titles = c("temperature", "dew point", "wind speed", "precipitation", "humidity", "pm2.5", "ozone")
for (i in 1:7) {
  plot(ny$t, ny[[cols[i]]], main=titles[i], xlab='time', pch=46)
}
plot(ny$t, ny$case, main='Infected cases', xlab='time', pch=46)

## Fit varying coefficient model
## Conditioner and predictors
u = covid.environment$t
cols = c('temp', 'wind', 'precipitation', 'humidity', 'pm25', 'ozone')
X = covid.environment[,cols]
X = cbind(1, X)
for (i in 2:dim(X)[2])
  X[,i] = (X[,i] - mean(X[,i])) / sd(X[,i])

# Fit varying coefficient model and visualize
y = covid.environment$log_case_avg
fit1 = vcm.asf.predictor.specific(X, y, u)
fit2 = vcm.asf.equidistant(X, y, u)
B1 = fit1$coef(u)
B2 = fit2$coef(u)
t = covid.environment$time
titles = c('intercept', 'temperature', 'wind speed', 'precipitation',
           'humidity', 'pm2.5', 'ozone')
par(mfrow=c(2,4), mar=c(4,2,2,1.5))
for (i in 1:7) {
  ind = order(u)
  lower = min(B1[,i], B2[,i]) - 0.2
  upper = max(B1[,i], B2[,i]) + 0.2
  plot(t[ind], B1[ind, i], lwd=2.0, type='n', xlab='time', ylab=paste0('B', i), ylim=c(lower, upper), main=titles[i])
  grid(6, NA, col='grey', lwd = 2)
  lines(t[ind], B1[ind, i], lwd=2)
  lines(t[ind], B2[ind, i], lty=2, lwd=2)
}

# Use rolling window to predict Covid cases
# fitting window size: at least 1 year, validation size: 1 week
tmin = as.Date("2021-02-28")
tmax = as.Date("2021-09-30")
validation = 7
start = tmin
end = start + validation
y = covid.environment$log_case_avg
ypred = covid.environment[covid.environment$time > tmin, 'log_case_avg']
f1 = c()
f2 = c()
count = 1
while (start <= tmax) {
  training = which(covid.environment$time <= start)
  testing = which((covid.environment$time > start) & (covid.environment$time <= end))
  Xtraining = X[training,]
  utraining = u[training]
  ytraining = y[training]
  Xtesting = X[testing,]
  utesting = u[testing]
  fit1 = vcm.asf.predictor.specific(Xtraining, ytraining, utraining, boundary = c(min(utraining), max(utesting)))
  fit2 = vcm.asf.equidistant(Xtraining, ytraining, utraining, boundary = c(min(utraining), max(utesting)))
  f1 = c(f1, fit1$predict(Xtesting, utesting))
  f2 = c(f2, fit2$predict(Xtesting, utesting))
  start = start + validation
  end = end + validation
  count = count + 1
}
# MSE for equidistant and predictor-specific spline fitting
print(c(mean((ypred - f2)^2), mean((ypred - f1)^2)))

# Fit varying coefficient model with different lags and compare RMSE
rmse = rep(0, 22)
y = covid.environment$log_case_avg
fit = vcm.asf.predictor.specific(X, y, u)
rmse[1] = sqrt(mean((y - fit$predict(X, u))^2))
for (lag in 1:21) {
  col = paste0('case_avg_next_', lag)
  covid.environment[which(covid.environment[[col]] < 1), col] = 1
  y = log(covid.environment[[col]])
  fit = vcm.asf.predictor.specific(X, y, u)
  rmse[lag+1] = sqrt(mean((y - fit$predict(X, u))^2))
}
plot(seq(0,21,1), rmse, xlab='lag', ylab='RMSE', main='RMSE for different lags', type='l')