# 0.a Import modules, global parameters
library(Matrix)

p = 6  # Dimension
n = 75 # Number of data points
alpha = 0.1  # Significance level
beta0 = 0  # "Unknown" intercept
beta = c(0.5, 0.25, 0.15, 0.6, 2.5, 1.2) # "Unknown" weights


# 0.b Code helper functions
LogLogitLikelihood <- function(b, y, x) {  # append a column of 1s for the x-values, to take into account intercept
  x = cbind(rep(1, n), x)
  out = 0
  for (i in 1:length(y)) {
    dot = x[i,] %*% b
    out = out + y[i] * dot - log(1 + exp(dot))
  }
  return (out)
}

LogitLikelihood <- function(b, y, x) {
  return (exp(LogLogitLikelihood(b, y, x)))
}

LRTfunc <- function(t, e_vect, intercept) {
  L1 = LogLogitLikelihood(estimator + t*e_vect, y, x)
  L0 = LogLogitLikelihood(estimator, y, x)
  return(intercept - 2 * (L0 - L1))
}

EllipsoidVolume <- function(ellipsoid, intercept) {
  V_sphere = pi^((p+1)/2) / gamma((p + 3)/2)
  return (sqrt(exp(sum(log(eigen(ellipsoid)$values))) * intercept^(p+1)) * V_sphere)
}


# 1. Generate the data
x = matrix(rnorm(n*p), n, p)
prob = 1/(1 + exp(-x%*%beta - beta0))
y = rbinom(n,1,prob)
model = glm(y ~x, family = binomial)


# 2. Create initial confidence region based on Wald's theorem
estimator = model$coefficients
ellipsoid = solve(vcov(model))
CritReg <- function(theta){  # theta has to be of dim p+1. Allows you to calculate p-values
  return((estimator - theta) %*% inv_finfo %*% (estimator - theta))
}


# 3. Optimise along e-vector axes using LRT. We start along a single axis only
# If qchisq(1 - alpha, p+1) - 2 * (L0 - L1) < 0, then do nothing.

Stretch <- function(direction, scale) {
  y = direction + c(-norm(direction, type="2"), rep(0, length(direction) - 1))
  w = y/norm(y, type="2")
  Q_w = diag(length(direction)) - 2 * w %*% t(w)
  S = diag(length(direction))
  S[1,1] = scale
  return (Q_w %*% S %*% Q_w)
}

OptimiseArea <- function(ellipsoid, y, x) {  # Returns the appropriate scaling matrix
  crit = qchisq(1 - alpha, p+1)
  res = diag(p+1)
  for (i in 1:(p+1)) {
    e_vect = eigen(ellipsoid)$vectors[i, ]
    x0 = (crit/(t(e_vect) %*% ellipsoid %*% e_vect))^0.5
    starting_v = estimator + x0*e_vect
    L1 = LogLogitLikelihood(starting_v, y, x)
    L0 = LogLogitLikelihood(estimator, y, x)
    if (is.nan(crit - 2 * (L0 - L1))){
      next
    } else if(crit - 2 * (L0 - L1) > 0){
      root = uniroot(LRTfunc, interval=c(x0, 100 + x0), e_vect=e_vect, intercept=crit)$root
      res = res %*% Stretch(e_vect, root/x0)
    } else {
      # I can't get the "shrink area if <0" part to work for some reason
      # root = uniroot(LRTfunc, interval=c(-100 + x0, x0 + 100), e_vect=e_vect, intercept=crit)$root
      # res = res %*% Stretch(e_vect, root/x0)
      
      next
    }
  }
  return(res)
}


# 4. Compare pre- and post- optimisation
JudgeOptim <- function(ellipsoid, y, x, estimator, b0, b, crit) {
  res = c(0, 0)  # First zero: preoptimisation makes true value in, second: optimisation makes true value in
  transform_ellipsoid = solve(OptimiseArea(ellipsoid, y, x))
  pre_factor = (estimator - c(beta0, beta))
  if ((t(pre_factor) %*% ellipsoid %*% pre_factor) < crit) {
    res[1] = 1
  } 
  post_factor = transform_ellipsoid %*% (estimator - c(beta0, beta))
  if ((t(post_factor) %*% ellipsoid %*% post_factor) < crit) {
    res[2] = 1
  } 
  post_vol = EllipsoidVolume(transform_ellipsoid %*% ellipsoid %*% transform_ellipsoid, crit)
  pre_vol = EllipsoidVolume(ellipsoid, crit)
  return(c(res, (post_vol - pre_vol)/pre_vol))
}


# 5. Now for the 10,000 iterations!
options(warn=-1)
crit = qchisq(1-alpha, p+1)
pre_optim = 0
post_optim = 0
n_iter = 10000
volume_change = c()
for (i in 1:n_iter) {
  # Progress Bar
  if (i %% 50 == 0) {
    gc()
    print(i)
  }

  # Generate data
  x = matrix(rnorm(n*p), n, p)
  prob = 1/(1 + exp(-x%*%beta - beta0))
  y = rbinom(n,1,prob)
  model = glm(y ~x, family = binomial)

  # Create initial confidence region
  estimator = model$coefficients
  ellipsoid = solve(vcov(model))

  # Now we check if optimisation makes it better!
  judge = JudgeOptim(ellipsoid, y, x, estimator, beta0, beta, crit)
  pre_optim = pre_optim + judge[1]
  post_optim = post_optim + judge[2]
  volume_change = c(volume_change, judge[3])
}

pre_optim
post_optim
n_iter
