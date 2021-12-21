# Application to COVID-19 case numbers and environmental factors dataset
library(vcmasf)

# Data preparation 
covid.environment[which(covid.environment$case_avg < 1), 'case_avg'] = 1
covid.environment['log_case_avg'] = log(covid.environment['case_avg'])
covid.environment['time'] = as.Date(covid.environment$time, format='%Y-%m-%d')
cols = c('temp', 'wind', 'precipitation', 'humidity', 'pm25', 'ozone')
covid.environment['year'] = unclass(as.POSIXlt(covid.environment$time))$year
covid.environment['year'] = covid.environment['year'] - 120
covid.environment['yday'] = unclass(as.POSIXlt(covid.environment$time))$yday 
covid.environment[which(covid.environment$year == 0), 'yday'] = covid.environment[which(covid.environment$year == 0), 'yday'] / 366
covid.environment[which(covid.environment$year != 0), 'yday'] = covid.environment[which(covid.environment$year != 0), 'yday'] / 365
covid.environment['t'] = covid.environment['year'] + covid.environment['yday']

# Conditioner and predictors
u = covid.environment$t
X = covid.environment[,cols]
X = cbind(1, X)
for (i in 2:dim(X)[2])
  X[,i] = (X[,i] - mean(X[,i])) / sd(X[,i])

# Fit varying coefficient model 
y = covid.environment$log_case_avg
fit = vcm.asf.twostep(X, y, u)

# Fit varying coefficient model with lag = 4
covid.environment[which(covid.environment$case_avg_next_4 < 1), 'case_avg_next_4'] = 1
covid.environment['log_case_avg_next_4'] = log(covid.environment['case_avg_next_4'])
y_lag = covid.environment$log_case_avg_next_4
fit_lag = vcm.asf.twostep(X, y_lag, u)

# Visualize coefficients
X = X[order(u),]
u = sort(u)
b = fit$coef(u)
b_lag = fit_lag$coef(u)
xticks = seq(as.Date("2020-06-01", format='%Y-%m-%d'), as.Date("2021-06-30", format='%Y-%m-%d'), "6 months")
xticksshow = c('2020-06', '2020-12', '2021-06')
titles = c('temperature', 'dew point', 'wind speed', 'precipitation', 'humidity', 'pm2.5', 'ozone', 'Infected cases')
par(mfrow=c(2, 4))
lower = min(cbind(b[,1], b_lag[,1]))
upper = max(cbind(b[,1], b_lag[,1]))
plot(covid.environment$time, b[,1], type='l', xlab='time', ylab='beta1', main='intercept', xaxt='n', ylim=c(lower, upper))
lines(covid.environment$time, b_lag[,1], lty=2, xlab='time', ylab=paste0('beta', i), main=titles[i-1], xaxt='n', )
axis(side=1, at=xticks, labels=xticksshow)
for (i in 2:dim(X)[2]) {
  lower = min(cbind(b[,i], b_lag[,i]))
  upper = max(cbind(b[,i], b_lag[,i]))
  plot(covid.environment$time, b[,i], type='l', xlab='time', ylab=paste0('beta', i), main=titles[i-1], xaxt='n', ylim=c(lower, upper))
  lines(covid.environment$time, b_lag[,i], lty=2, xlab='time', ylab=paste0('beta', i), main=titles[i-1], xaxt='n')
  axis(side=1, at=xticks, labels=xticksshow)
}