# setwd("~/Documents/phd_Qua/qua4/experiments/simulation")
source("../inference/util.r")
source("../inference/MM_estimator.r")
source("../data/synthetic_dat.r")
# setup julia 
library(JuliaCall)
julia<- julia_setup()
# source esl_lasso implemented in julia 
julia_call("include", "../inference/esl_lasso.jl") 




# by default only repeat step1-3 once
ESL_LASSO <-function(dat, beta_true  = c(1,1.5, 2,1, 0,0,0,0)){
  d = dim(dat)[2]
  X = as.matrix(dat[,2:d])
  
  # get MM_estimator for initial value
  MM_fit = MM_estimator(dat, beta_true)
  init_est = c(MM_fit$coef)
  
  # call julia function to run esl_lasso
  dat = as.matrix(dat)
  julia_assign("beta0", init_est)
  julia_assign("dat",dat)
  coef = julia_eval("esl_lasso(beta0, dat)")  
  
  # measures of performance
  psr = PSR(as.vector(coef),beta_true)
  nsr = NSR(as.vector(coef),beta_true)
  me = ME(X, as.vector(coef),beta_true)
  
  result = list(coef = coef, PSR = psr, NSR = nsr, ME = me)
  return(result)
}





# ESL_LASSO with new gamma selection step
our_ESL_LASSO <-function(dat, beta_true  = c(1,1.5, 2,1, 0,0,0,0)){
  d = dim(dat)[2]
  X = as.matrix(dat[,2:d])
  
  # get MM_estimator for initial value
  MM_fit = MM_estimator(dat, beta_true)
  init_est = c(MM_fit$coef)
  
  # call julia function to run esl_lasso
  dat = as.matrix(dat)
  julia_assign("beta0", init_est)
  julia_assign("dat",dat)
  coef = julia_eval("our_esl_lasso(beta0, dat)")  
  
  # measures of performance
  psr = PSR(as.vector(coef),beta_true)
  nsr = NSR(as.vector(coef),beta_true)
  me = ME(X, as.vector(coef),beta_true)
  
  result = list(coef = coef, PSR = psr, NSR = nsr, ME = me)
  return(result)
}





# dat = df1_gen(800)
# our_ESL_LASSO(dat)



