group.lasso.prepare = function(X, y, u, w, Knots, nums, boundary, degree) {
  my = mean(y)
  sdy = mean((y - my)^2)
  y = (y - my) / sdy
  if (is.null(w))
    w = rep(1, dim(X)[2])
  else
    w = w / mean(w)
  if (is.null(boundary))
    boundary = range(u)
  p = dim(X)[2]
  n = dim(X)[1]
  rmsx = rep(0, p)
  for (i in 1:p) {
    rmsx[i] = mean(X[,i]^2)
    X[,i] = X[,i] / rmsx[i]
  }
  
  if (is.null(dim(Knots))) {
    nums = rep(length(Knots), p)
    Knots = t(matrix(rep(Knots, p), ncol=p))
  }
  
  U = matrix(nrow=n, ncol=0)
  start_ind = c()
  end_ind = c()
  ind = 1
  for (i in 1:p) {
    if (nums[i] > 0)
      tmp = bSpline(u, knots=Knots[i, 1:nums[i]], degree=degree, Boundary.knots=boundary, intercept=TRUE)
    else
      tmp = bSpline(u, knots=NULL, degree=degree, Boundary.knots=boundary, intercept=TRUE)
    tmp = gramSchmidt(X[,i] * gramSchmidt(tmp)$Q)$Q * sqrt(n)
    U = cbind(U, tmp) 
    start_ind = c(start_ind, ind)
    end_ind = c(end_ind, ind + nums[i] + degree)
    ind = ind + nums[i] + degree + 1
  }
  
  return (list(y=y, U=U, start_ind=start_ind, end_ind=end_ind, nums=nums, degree=degree, w=w, Knots=Knots))
}

group.lasso.iteration = function(data, lambda, tol, max_iter) {
  y = data$y
  n = length(y)
  U = data$U
  start_ind = data$start_ind
  end_ind = data$end_ind
  nums = data$nums
  degree = data$degree
  w = data$w
  p = length(nums)
  gamma = matrix(0, nrow=p, ncol=max(nums) + degree + 1)
  iter = 0
  diff = 2 * tol
  r = y
  
  while ((iter < max_iter) & (diff > tol)) {
    diff = 0
    iter = iter + 1
    for (i in 1:p) {
      gamma_old = gamma[i, 1:(nums[i] + degree + 1)]
      z = colMeans(r * U[,start_ind[i]:end_ind[i]]) + gamma_old 
      tmp = 1 - lambda * w[i] / sqrt(sum(z^2))
      tmp[tmp < 0] = 0
      gamma_new = tmp * z
      r = r - as.numeric(U[,start_ind[i]:end_ind[i]] %*% (gamma_new  - gamma_old))
      gamma[i, 1:(nums[i] + degree + 1)] = gamma_new
      diff = diff + sum((gamma_new - gamma_old)^2)
    }
    diff = sqrt(diff)
  }
  
  group_size = rowMeans(gamma^2)
  df = 0
  f = rep(0, n)
  for (i in 1:p)
  {
    if (group_size[i] > 0) {
      f = f + as.numeric(U[,start_ind[i]:end_ind[i]] %*% gamma[i, 1:(nums[i] + degree + 1)])
      df = df + nums[i] + degree + 1
    }
  }
  
  bic = log(sum((f - y)^2)) + log(n)*df/n
  select_inds = c()
  for (i in 1:p) {
    if (group_size[i] > 0)
      select_inds = c(select_inds, i)
  }
  return (list(select_inds=select_inds, bic=bic, gamma=gamma))
}

group.lasso = function(X, y, u, w, Knots, nums, boundary=NULL, degree=3, lambda=1.0, tol=1e-8, max_iter=20) {
  data = group.lasso.prepare(X, y, u, w, Knots, nums, boundary, degree)
  return (group.lasso.iteration(data, lambda, tol, max_iter))
}


group.lasso.bic = function(X, y, u, w, Knots, nums, boundary=NULL, degree=3, total_lambdas=100, tol=1e-8, max_iter=20) {
  data = group.lasso.prepare(X, y, u, w, Knots, nums, boundary, degree)
  y = data$y
  U = data$U
  start_ind = data$start_ind
  end_ind = data$end_ind
  degree = data$degree
  nums = data$nums
  w = data$w
  lam_max = 0
  p = dim(X)[2]
  for (i in 1:p) 
    lam_max = max(lam_max, sqrt(mean((y * U[, start_ind[i]:end_ind[i]])^2) / w[i]))
 
  lambdas = lam_max * exp(3*log(10) * c(0:(total_lambdas-1))/(total_lambdas-1)-3*log(10))
  lambda_star = 0
  bic = Inf 
  select_inds = c()
  gamma = matrix(0, nrow=p, ncol=max(nums)+degree+1)
  for (i in 1:total_lambdas) {
    tmp = group.lasso.iteration(data, lambdas[i], tol, max_iter)
    if (tmp$bic < bic) {
      bic = tmp$bic
      select_inds = tmp$select_inds
      lambda_star = lambdas[i]
      gamma = tmp$gamma
    }
  }
  
  if (is.null(select_inds))
    return (list(select_inds=select_inds, gamma=NULL, nums=NULL, Knots=NULL))
  
  return(list(select_inds=select_inds, gamma=matrix(gamma[select_inds,], nrow=length(select_inds)), 
              nums=nums[select_inds], Knots=data$Knots[select_inds,]))
}

adaptive.group.lasso.bic = function(X, y, u, Knots, nums, boundary=NULL, degree=3, total_lambdas=100, tol=1e-8, max_iter=20) {
  glasso = group.lasso.bic(X, y, u, NULL, Knots, nums, boundary, degree, total_lambdas, tol, max_iter)
  if (is.null(glasso$select_inds))
    return (glasso$select_inds)
  
  inds = glasso$select_inds
  Knots = glasso$Knots
  nums = glasso$nums
  num_selected = length(inds)
  n = dim(X)[1]
  X_new = matrix(0, nrow=n, ncol=num_selected)
  w_new = rep(0, num_selected)
  for (i in 1:num_selected) {
    ind = glasso$select_inds[i]
    X_new[,i] = X[,ind]
    w_new[i] = 1 / sqrt(sum(glasso$gamma[i,1:(glasso$nums[i]+degree+1)]^2))
  }
  aglasso = group.lasso.bic(X_new, y, u, w_new, Knots, nums, boundary, degree, total_lambdas, tol, max_iter)
  aglasso$select_inds = inds[aglasso$select_inds]
  
  return (aglasso$select_inds)
}
