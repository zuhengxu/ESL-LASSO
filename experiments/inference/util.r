
# positive selection rate
PSR<- function(beta, beta_true ){
  ind1 = which(beta_true !=0)
  ps = beta[ind1]
  psr = length(ps[ps!=0])/length(ind1)
  return(psr)
}

# non-causal selection rate
NSR <-function(beta, beta_true){
  ind0 = which(beta_true == 0)
  ns  = beta[ind0]
  nsr = length(ns[ns==0])/length(ind0)
  return(nsr)
}

# model error (beta- beta0)^T E[xx^T] (beta- beta0)
ME <-function(X, beta, beta_true){
  X  =as.matrix(X)
  E = t(X)%*%X/dim(X)[1]
  diff = as.matrix(beta - beta_true)
  me = t(diff)%*% E %*%diff
  return(me)
}

