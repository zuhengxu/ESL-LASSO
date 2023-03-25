library(pense)
# setwd("~/Documents/phd_Qua/qua4/experiments/inference")
source("../inference/util.r")

# input data with Y = dat[,1], X = dat[,2:]
pense_lasso <-function(dat, cv_k = 4, beta_true =  c(1, 1.5, 2, 1, 0,0,0,0)){
  d = dim(dat)[2]
  X = as.matrix(dat[,2:d])
  
  # CV: choose lambda that leads to minimal mean square prediction error
  # don't include intercept
  cv_results <- pense_cv(X, dat$Y, alpha = 1, cv_k = cv_k, intercept = FALSE)
  coef = t(as.matrix(coef(cv_results, lambda = 'min')[2:9]))
  
  # measures of performance
  psr = PSR(as.vector(coef),beta_true)
  nsr = NSR(as.vector(coef),beta_true)
  me = ME(X, as.vector(coef),beta_true)
  
  result = list(coef = coef, PSR = psr, NSR = nsr, ME = me)
  return(result)
}


# #Ex
# dat = df1_gen(100)
# pense_fit = pense_lasso(dat)
