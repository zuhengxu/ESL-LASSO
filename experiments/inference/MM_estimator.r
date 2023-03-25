library(robustbase)
# setwd("~/Documents/phd_Qua/qua4/experiments")
source("../inference/util.r")


MM_estimator<-function(dat, beta_true =  c(1, 1.5, 2, 1, 0,0,0,0)){
  d = dim(dat)[2]
  X = as.matrix(dat[,2:d])
  
  # no intercept
  fit = lmrob(Y ~ .-1, data = dat, method = "MM") 
  #coefficients
  coef = t(as.matrix(fit$coefficients))
  # measures of performance
  psr = PSR(as.vector(coef),beta_true)
  nsr = NSR(as.vector(coef),beta_true)
  me = ME(X, as.vector(coef),beta_true)
  
  result = list(coef = coef, PSR = psr, NSR = nsr, ME = me)
  return(result)
}


 



