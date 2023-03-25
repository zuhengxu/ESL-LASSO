# setwd("~/Documents/phd_Qua/qua4/experiments")
source("../data/synthetic_dat.r")
source("../inference/PENSE_LASSO.r")
source("../inference/ESL_LASSO.r")
# julia setup 
library(JuliaCall)
julia<- julia_setup()
# source julia file (esl_lasso implementation)
julia_call("include", "../inference/esl_lasso.jl") 




args = commandArgs(trailingOnly =TRUE)
N =  as.integer(args[1]) # sample size
ID = as.integer(args[2]) # ID for the trials
f = get(args[3]) # dataset gen function
set.seed(ID)


Names = c("method",'PSR', 'NSR', 'ME' , 'x1', 'x2', 'x3', 'x4', 'x5', 'x6','x7','x8')
df_pense = matrix(, nrow = 20, ncol = 12)
df_esl = matrix(, nrow = 20, ncol = 12)
df_our_esl = matrix(, nrow = 20, ncol = 12)

for (i in 1:20){
  dat = f(N)
  pense_fit = pense_lasso(dat)
  MM_fit = MM_estimator(dat)
  esl_fit = ESL_LASSO(dat) 
  our_esl_fit = our_ESL_LASSO(dat) 
  
  df_pense[i,] = c("PENSE-LASSO",pense_fit$PSR, pense_fit$NSR, pense_fit$ME, pense_fit$coef)
  df_esl[i,] = c("ESL-LASSO", esl_fit$PSR, esl_fit$NSR, esl_fit$ME, esl_fit$coef)
  df_our_esl[i,] = c("ESL-LASSO(ours)", our_esl_fit$PSR, our_esl_fit$NSR, our_esl_fit$ME, our_esl_fit$coef)
}

#append all result into 
df_result = as.data.frame(rbind(df_pense,df_esl, df_our_esl))
colnames(df_result) = Names

#write results table
write.table(df_result, paste0("../results/",N , "_",args[3],".csv"), append = TRUE,
            row.names = FALSE, col.names = !file.exists(paste0("../results/", N, "_",args[3],".csv")))
