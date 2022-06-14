## Update the code so the repeat u values are in the same slice.
knots.selection = function(X, y, u, ms = ceiling(sqrt(length(y))), lambda0 = 3.0, Knots=c()) {
  n = dim(X)[1]
  p = dim(X)[2]
  y = (y - mean(y)) / sd(y)
  X = t(t(X) / sqrt(colMeans(X^2)))
  
  y = y[order(u)]
  X = as.matrix(X[order(u),])
  u = u[order(u)]
  if (length(Knots) == 0)
    knots_index = rev(knots_selection_cpp(X, y, ms, lambda0, u, u))
  else
    knots_index = rev(knots_selection_cpp(X, y, 0, lambda0, Knots, u))
  s = length(knots_index) - 1
  if (s == 0) 
    knots = NULL
  else {
    knots = rep(0, s)
    if (length(Knots) == 0) {
      for (i in 1:s)
        knots[i] = (u[knots_index[i]] + u[knots_index[i]+1]) / 2
    } else {
      for (i in 1:s) {
        for (j in 1:length(Knots)) {
          if (u[knots_index[i]] <= Knots[j]) {
            knots[i] = Knots[j]
            break
          }
        }
      }
    }
  }
  return (knots)
}

vcm.asf.onestep <- function(X, y, u, ms = ceiling(sqrt(length(y))), lambda0s=seq(0.1, 10, 0.1), boundary=NULL, degree=3, linear=TRUE, Knots=c())
{
  if (is.null(boundary)) 
    boundary = range(u)
  if (linear)
    X0 = cbind(X, X*u)
  else
    X0 = matrix(X, nrow=length(y))
  k = length(lambda0s)
  p = dim(X)[2]
  n = length(y)
  
  bic_min = Inf
  knots = NULL
  fit = NULL
  
  for (i in 1:k) {
    tmp = knots.selection(X0, y, u, ms, lambda0s[i], Knots)
    uspline = bSpline(u, knots=tmp, degree=degree, Boundary.knots=boundary, intercept=TRUE)
    bspline = matrix(0, nrow=n, ncol=0)
    for (j in 1:p) 
      bspline = cbind(bspline, X[,j] * uspline)
    tmp_fit = lm(y ~ bspline-1)
    bic = log(mean((y - tmp_fit$fitted.values) ** 2)) + dim(bspline)[2]*log(n) / n
    if (bic < bic_min) {
      bic_min = bic
      lambda0 = lambda0s[i]
      knots = tmp
      fit = tmp_fit
    }
  }

  predict <- function(Xnew, unew) {
    uspline_new = bSpline(unew, knots=knots, degree=degree, Boundary.knots=boundary, intercept=TRUE)
    bspline_new = matrix(0, nrow=length(unew), ncol=0)
    for (j in 1:p)
      bspline_new = cbind(bspline_new, Xnew[,j] * uspline_new)
    return (c(bspline_new %*% fit$coefficients))
  }
  
  coef <- function(unew) {
    coefs = matrix(0, nrow=length(unew), ncol=p)
    uspline_new = bSpline(unew, knots=knots, degree=degree, Boundary.knots=boundary, intercept=TRUE)
    ps = dim(uspline_new)[2]
    for (j in 1:p)
      coefs[,j] = c(uspline_new %*% fit$coefficients[(1+(j-1)*ps):(j*ps)])
    
    return (coefs)
  }
  
  return(list(knots=knots, predict=predict, coef=coef))
}

vcm.asf.twostep = function(X, y, u, ms = ceiling(sqrt(length(y))), lambda0s=seq(0.1, 10, 0.1), boundary=NULL, degree=3, max_iter=100, linear=TRUE, Knots=c()) {
  Knots_orig = Knots
  res = vcm.asf.onestep(X, y, u, ms, lambda0s, boundary, degree, linear=linear, Knots_orig)
  coef = res$coef(u)
  if (is.null(boundary)) 
    boundary = range(u)
  n = dim(X)[1]
  p = dim(X)[2]
  Knots = matrix(0, nrow=p, ncol=ceiling(n/ms))

  Splines = list()
  uspline = bSpline(u, knots=res$knots, degree=degree, Boundary.knots=boundary, intercept=TRUE)
  bspline = matrix(0, nrow=n, ncol=0)
  ind = c(0)
  nums = rep(0, p)
  for (i in 1:p) {
    nums[i] = length(res$knots)
    if (nums[i] > 0)
      Knots[i,1:nums[i]] = res$knots
    Splines[[i]] = uspline * X[,i]
    ind = c(ind, dim(Splines[[i]])[2])
    bspline = cbind(bspline, Splines[[i]])
  }
  ind = cumsum(ind)
  fit = matrix(0, nrow=n, ncol=p)
  fit_lm = lm(y~bspline-1)
  ## Remove na coefficients
  fit_lm$coefficients[is.na(fit_lm$coefficients)] = 0
  for (i in 1:p) {
    fit[,i] = c(Splines[[i]] %*% fit_lm$coefficients[(ind[i]+1):(ind[i+1])])
  }
  min_bic = log(mean((y - fit_lm$fitted.values) ** 2)) + dim(bspline)[2] * log(n)/n 
  iter = 0
  while (iter < max_iter) {
    iter = iter + 1
    min_feature = 0
    min_knots = c()
    for (j in 1:p) {
      residual = y - rowSums(fit) + fit[,j]
      prev = c()
      for (lambda0 in lambda0s)
      {
        tmp = vcm.asf.onestep(matrix(X[,j], ncol=1), residual, u, ms, c(lambda0), boundary, degree, linear=linear, Knots=Knots_orig) 
        if (length(prev) == length(tmp$knots)) {
          if (sum((prev - tmp$knots)^2) == 0)
            next
        }
        prev = tmp$knots
        base_spline = matrix(0, nrow=n, ncol=0)
        for (k in 1:p) {
          if (k == j)
            base_spline = cbind(base_spline, bSpline(u, knots=tmp$knots, degree=degree, Boundary.knots=boundary, intercept=TRUE) * X[,j])
          else
            base_spline = cbind(base_spline, Splines[[k]])
        }
        tmp_lm = lm(y~base_spline-1)
        bic = log(mean((y - tmp_lm$fitted.values) ** 2)) + dim(base_spline)[2] * log(n)/n
        if (bic < min_bic) {
          min_bic = bic
          min_feature = j
          min_knots = tmp$knots
        }
      }
    }
    if (min_feature == 0)
      break
    nums[min_feature] = length(min_knots)
    Knots[min_feature,] = 0
    if (nums[min_feature] > 0)
      Knots[min_feature, 1:nums[min_feature]] = min_knots
    Splines[[min_feature]] = bSpline(u, knots=min_knots, degree=degree, Boundary.knots=boundary, intercept=TRUE) * X[,min_feature]
    bspline = matrix(0, nrow=n, ncol=0)
    ind = c(0)
    for (i in 1:p) {
      ind = c(ind, dim(Splines[[i]])[2])
      bspline = cbind(bspline, Splines[[i]])
    }
    ind = cumsum(ind)
    fit = matrix(0, nrow=n, ncol=p)
    fit_lm = lm(y~bspline-1)
    ## Remove na coefficients
    fit_lm$coefficients[is.na(fit_lm$coefficients)] = 0
    for (i in 1:p) {
      fit[,i] = c(Splines[[i]] %*% fit_lm$coefficients[(ind[i]+1):(ind[i+1])])
    }
  }
  knots_list = vector(mode='list', length=p)
  for (i in 1:p) {
    if (nums[i] > 0)
      knots_list[[i]] = Knots[i, 1:nums[i]]
  }
  
  predict <- function(Xnew, unew) {
    N = dim(Xnew)[1]
    Xnewspline = matrix(0, nrow=N, ncol=0)
    
    for (i in 1:p) {
      start = ind[i]+1
      end = ind[i+1]
      if (nums[i] > 0)
        Xnewspline = cbind(Xnewspline, bSpline(unew, knots=Knots[i, 1:nums[i]], degree=degree, Boundary.knots=boundary, intercept=TRUE) * Xnew[,i])
      else
        Xnewspline = cbind(Xnewspline, bSpline(unew, knots=NULL, degree=degree, Boundary.knots=boundary, intercept=TRUE) * Xnew[,i])
    }
    return (c(Xnewspline %*% fit_lm$coefficients))
  }
  
  coef <- function(unew) {
    N = length(unew)
    b = matrix(0, ncol=p, nrow=N)
    start = 1
    for (i in 1:p) {
      if (nums[i] > 0)
        tmp = bSpline(unew, knots=Knots[i, 1:nums[i]], degree=degree, Boundary.knots=boundary, intercept=TRUE)
      else
        tmp = bSpline(unew, knots=NULL, degree=degree, Boundary.knots=boundary, intercept=TRUE)
      b[,i] = tmp %*% fit_lm$coefficients[start:(start+dim(tmp)[2]-1)]
      start = start+dim(tmp)[2]
    }
    return (b)
  }
  
  return(list(knots=knots_list, predict=predict, coef=coef))
}

vcm.asf.equidistant <- function(X, y, u, nknots=seq(1, 9), boundary=NULL, degree=3)
{
  if (is.null(boundary)) 
    boundary = range(u)
  k = length(nknots)
  p = dim(X)[2]
  n = length(y)
  
  bic_min = Inf
  knots = NULL
  fit = NULL
  
  for (i in 1:k) {
    qs = seq(1, i) / (i+1)
    tmp = quantile(u, qs)
    uspline = bSpline(u, knots=tmp, degree=degree, Boundary.knots=boundary, intercept=TRUE)
    bspline = matrix(0, nrow=n, ncol=0)
    for (j in 1:p) 
      bspline = cbind(bspline, X[,j] * uspline)
    tmp_fit = lm(y ~ bspline-1)
    bic = log(mean((y - tmp_fit$fitted.values) ** 2)) + dim(bspline)[2]*log(n) / n
    if (bic < bic_min) {
      bic_min = bic
      knots = tmp
      fit = tmp_fit
    }
  }
  
  predict <- function(Xnew, unew) {
    uspline_new = bSpline(unew, knots=knots, degree=degree, Boundary.knots=boundary, intercept=TRUE)
    bspline_new = matrix(0, nrow=length(unew), ncol=0)
    for (j in 1:p)
      bspline_new = cbind(bspline_new, Xnew[,j] * uspline_new)
    return (c(bspline_new %*% fit$coefficients))
  }
  
  coef <- function(unew) {
    coefs = matrix(0, nrow=length(unew), ncol=p)
    uspline_new = bSpline(unew, knots=knots, degree=degree, Boundary.knots=boundary, intercept=TRUE)
    ps = dim(uspline_new)[2]
    for (j in 1:p)
      coefs[,j] = c(uspline_new %*% fit$coefficients[(1+(j-1)*ps):(j*ps)])
    
    return (coefs)
  }
  
  return(list(knots=knots, predict=predict, coef=coef))
}


knots.selection.marginal = function(X, y, u, ms=ceiling(sqrt(length(y))), lambda0s=seq(0.1, 10, 0.1), 
                                    boundary=NULL, degree=3, simple=TRUE, linear=TRUE, Knots=c()) {
  Knots_orig = Knots
  p = dim(X)[2]
  n = dim(X)[1]
  if (is.null(boundary))
    boundary = range(u)
  Knots = matrix(0, nrow=p, ncol=ceiling(n/ms))
  nums = rep(0, p)
  for (i in 1:p) 
  {
    tmp = vcm.asf.onestep(matrix(X[,i], ncol=1), y, u, ms, lambda0s, boundary, degree, linear=linear, Knots=Knots_orig)
    nums[i] = length(tmp$knots)
    if (nums[i] > 0)
      Knots[i, 1:nums[i]] = tmp$knots
  }
  return(list(Knots=Knots, nums=nums))
}

vcm.variable.selection = function(X, y, u, ms=ceiling(sqrt(length(y))), lambda0s=seq(0.1, 10, 0.1), boundary=NULL, degree=3, 
                                  total_glasso_penalties=100, tol=1e-8, max_iter=20, lam_min_ratio=1e-3, linear=TRUE, Knots=c())
{
  n = dim(X)[1]
  p = dim(X)[2]
  inds = seq(1, p)
  res = knots.selection.marginal(X, y, u, ms, lambda0s, boundary, degree, linear=linear, Knots=Knots)
  Knots = res$Knots
  nums = res$nums
  res = adaptive.group.lasso.bic(X, y, u, Knots, nums, boundary, degree, total_glasso_penalties, tol, max_iter, lam_min_ratio)
  glasso = res$glasso
  aglasso = res$aglasso
  if (length(glasso$select_inds) > 0)
    glasso$Knots = matrix(glasso$Knots, nrow=length(glasso$select_inds))
  if (length(aglasso$select_inds) > 0)
    aglasso$Knots = matrix(aglasso$Knots, nrow=length(aglasso$select_inds))
  
  predict <- function(Xnew, unew) {
    if (length(aglasso$select_inds) == 0)
      return (rep(mean(y), length(unew)))
    N = dim(Xnew)[1]
    Xnewspline = matrix(0, nrow=N, ncol=0)
    for (i in 1:length(aglasso$select_inds)) {
      if (aglasso$nums[i] > 0)
        Xnewspline = cbind(Xnewspline, bSpline(unew, knots=aglasso$Knots[i, 1:aglasso$nums[i]], degree=degree, Boundary.knots=boundary, 
                                               intercept=aglasso$intercepts[i]) * Xnew[,aglasso$select_inds[i]])
      else
        Xnewspline = cbind(Xnewspline, bSpline(unew, knots=NULL, degree=degree, Boundary.knots=boundary, 
                                               intercept=aglasso$intercepts[i]) * Xnew[,aglasso$select_inds[i]])
    }
    return (c(Xnewspline %*% aglasso$beta) + mean(y) - aglasso$alpha)
  }
  
  coef <- function(unew) {
    if (length(aglasso$select_inds) == 0)
      return (NULL)
    coefs = matrix(0, nrow=length(unew), ncol=length(aglasso$select_inds))
    start = 1
    for (i in 1:length(aglasso$select_inds)) {
      if (nums[i] > 0)
        tmp = bSpline(unew, knots=aglasso$Knots[i, 1:aglasso$nums[i]], degree=degree, Boundary.knots=boundary, 
                      intercept=aglasso$intercepts[i])
      else
        tmp = bSpline(unew, knots=NULL, degree=degree, Boundary.knots=boundary, 
                      intercept=aglasso$intercepts[i])
      coefs[,i] = tmp %*% aglasso$beta[start:(start+dim(tmp)[2]-1)]
      start = start+dim(tmp)[2]
    }
    return (coefs)
  }

  return (list(ind=aglasso$select_inds, predict=predict, coef=coef))
}

