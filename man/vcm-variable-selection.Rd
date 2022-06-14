\name{vcm.variable.selection}
\alias{vcm.variable.selection}
\title{Variable selection for sparse high dimension problem}
\description{
  This function performs variable selection for high dimension varying coefficient model.
}
\usage{
vcm.variable.selection(X, y, u, ms=ceiling(sqrt(length(y))), 
		       lambda0s=seq(0.1, 10, 0.1), boundary=NULL, 
		       degree=3, total_glasso_penalties=100, 
		       tol=1e-8, max_iter=20, lam_min_ratio=1e-3,
               linear=TRUE, Knots=c())
}
\arguments{
  \item{X}{a numeric matrix for predictors}
  \item{y}{a numeric vector for response}
  \item{u}{a numeric vector for conditioner}
  \item{ms}{minimum size of slices}
  \item{lambda0s}{a numeric vector for penalties}
  \item{boundary}{bounds for conditioner u to construct splines}
  \item{degree}{degree of splines}
  \item{total_glasso_penalties}{total number of penalties for adaptive group lasso variable selection}
  \item{tol}{maximum tolerance for adaptive group lasso update}
  \item{max_iter}{maximum iteration for adaptive group lasso update}
  \item{lam_min_ratio}{The ratio between the minimum and maximum penalties for adaptive group lasso}
  \item{linear}{If TRUE, use linear approximation within each segments; otherwise use constant approximation}
  \item{Knots}{If Knots is not c(), use the specify knots for fitting}
}
\value{
 The function returns a list with the following components
  \item{ind}{the selected variables}
  \item{coef}{a function that returns the varying coefficient given conditioner}
  \item{predict}{a function that returns the fitted values given predictors and conditioner} 
}
\references{
  See paper "Varying Coefficient Model via Adaptive Spline Fitting" for more details.
}
\examples{
  n = 500
  p = 10
  X = matrix(rnorm(n * p), ncol = p)
  u = runif(n)
  B0 = matrix(0, nrow=n, ncol=p)
  B0[,1] = sin(u * 2 * pi)
  B0[,2] = cos(u * 2 * pi)
  y = rowSums(X * B0) + rnorm(n)
  ind = vcm.variable.selection(X, y, u)
}
\author{Xufei Wang}

