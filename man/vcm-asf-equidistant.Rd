\name{vcm.asf.equidistant}
\alias{vcm.asf.equidistant}
\title{Spline fitting with equidistant knots}
\description{
  This function estimates the varying coefficient by spline fitting with equidistant knots.
}
\usage{
vcm.asf.equidistant(X, y, u, nknots=seq(1, 9), boundary=NULL, degree=3)
}
\arguments{
  \item{X}{a numeric matrix for predictors}
  \item{y}{a numeric vector for response}
  \item{u}{a numeric vector for conditioner}
  \item{nknots}{a numeric vector for number of knots}
  \item{boundary}{bounds for conditioner u to construct splines}
  \item{degree}{degree of splines}
}
\value{
  The function returns a list with the following components
  \item{knots}{a numeric vector of selected knots}
  \item{coef}{a function that returns the varying coefficient given conditioner}
  \item{predict}{a function that returns the fitted values given predictors and conditioner}
  
}
\references{
  See paper \href{https://arxiv.org/abs/2201.10063}{Varying Coefficient Model via Adaptive Spline Fitting} for more details.
}
\examples{
  n = 500
  p = 2
  X = matrix(rnorm(n * p), ncol = p)
  u = runif(n)
  B0 = matrix(0, nrow=n, ncol=p)
  B0[,1] = sin(u * 2 * pi)
  B0[,2] = cos(u * 2 * pi)
  y = rowSums(X * B0) + rnorm(n)
  fit = vcm.asf.equidistant(X, y, u)
  B = fit$coef(u)
  f = fit$predict(X, u)
}
\author{Xufei Wang}

