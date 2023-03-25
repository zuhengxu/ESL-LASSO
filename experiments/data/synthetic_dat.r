library(mvtnorm)
library(purrr)

# sampler from 2component gaussian mixture
rnormmix <- function(N, p1, mu1, mu2, omega1, omega2){
  r1 = rmvnorm(N, mu1, omega1)
  r2 = rmvnorm(N, mu2, omega2)
  Ind1 = rbernoulli(N, p1)
  Ind2  =1- Ind1
  return(Ind1*r1 + Ind2*r2) 
}


#1. influential points in the predictors
df1_gen <-function(N){
  #parameters
  beta_true = c(1, 1.5, 2, 1, 0,0,0,0)
  d = length(beta_true)
  I  =  matrix(1:d, d,d)
  J = 1:d
  omega2 = 0.5**abs(t(t(I) -J))
  omega1 = diag(d)
  mu1 =  rep(0, d)
  mu2 = rep(3, d)
  
  #datagen
  X1 = rnormmix(N, 0.8, mu1, mu2, omega1, omega2)
  Y1 = X1%*% matrix(beta_true) + rnorm(N)
  dat1 = cbind(Y1, X1)
  df1 = as.data.frame(dat1)
  colnames(df1) <- c('Y', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6','x7','x8')
  return(df1)
}


#2. influential points in  the response
df2_gen<-function(N){
  #parameters
  beta_true = c(1, 1.5, 2, 1, 0,0,0,0)
  d = length(beta_true)
  I  =  matrix(1:d, d,d)
  J = 1:d
  omega2 = 0.5**abs(t(t(I) -J))
  omega1 = diag(d)
  mu1 =  rep(0, d)
  mu2 = rep(3, d)
  
  #datagen
  X2 = rmvnorm(N, mu1, omega2)
  err2 =  rnormmix(N, 0.8, c(0), c(10), matrix(1, 1,1), matrix(36,1,1))
  Y2 = X2%*% matrix(beta_true) +err2
  dat2 = cbind(Y2, X2)
  df2 = as.data.frame(dat2)
  colnames(df2)<- c('Y', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6','x7','x8')
  return(df2)
}


#3. influential points in both the predictors and response
df3_gen <-function(N){
  #parameters
  beta_true = c(1, 1.5, 2, 1, 0,0,0,0)
  d = length(beta_true)
  I  =  matrix(1:d, d,d)
  J = 1:d
  omega2 = 0.5**abs(t(t(I) -J))
  omega1 = diag(d)
  mu1 =  rep(0, d)
  mu2 = rep(3, d)
  
  #datagen
  X3 = rnormmix(N, 0.8, mu1, mu2, omega1, omega2)
  err3 = rcauchy(N)
  Y3 = X3%*%matrix(beta_true) + err3
  dat3 = cbind(Y3, X3)
  df3 = as.data.frame(dat3)
  colnames(df3)<- c('Y', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6','x7','x8')
  return(df3)
}


#4. test dataset 
df_test_gen<- function(N){
  #parameters
  beta_true = c(1, 1.5, 2, 1, 0,0,0,0)
  d = length(beta_true)
  I  =  matrix(1:d, d,d)
  J = 1:d
  omega2 = 0.5**abs(t(t(I) -J))
  mu1 =  rep(0, d)
  
  #datagen
  X = rmvnorm(N, mu1, omega2)
  err = rcauchy(N)
  Y = X%*%matrix(beta_true) + err
  dat = cbind(Y, X)
  df = as.data.frame(dat)
  colnames(df)<- c('Y', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6','x7','x8')
  return(df)
}

