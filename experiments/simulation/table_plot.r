setwd("~/Documents/phd_Qua/qua4/experiments")
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(wesanderson)


# read the simulation result and show violin plot
N  = 400
ID = 2
df  = read_delim(paste0("results/",N ,"_", ID, "_gen.csv"), escape_double  =FALSE, delim = " ")
df = as.tibble(df)
df %>% select('method', 'PSR', 'NSR', 'ME') %>% 
  filter(method %in% c("ESL-LASSO(ours)", "PENSE-LASSO") )%>% 
  ggplot(aes(x = method , y = ME, fill = method )) + geom_violin(trim = FALSE)+
  geom_boxplot(width = .1)+ 
  theme(legend.position="none", text = element_text(size=20)) +
  ggtitle(paste0("N = ", N, " , Simulation Setting = ", ID))+
  scale_fill_brewer(palette="BuPu")



# ME violin plot
p<- list(); i = 0
for(N in c(100,400, 800, 1200)){
  for (ID in c(1,2,3)){
    i = i+1
    df  = read_delim(paste0("results/",N ,"_df", ID, "_gen.csv"), escape_double  =FALSE, delim = " ")
    df = as.tibble(df)
    p[[i]] =  df %>% select('method', 'PSR', 'NSR', 'ME') %>% 
      filter(method %in% c("ESL-LASSO(ours)", "PENSE-LASSO") )%>% 
      ggplot(aes(x = method , y = ME, fill = method )) + geom_violin(trim = FALSE)+
      geom_boxplot(width = .1)+ 
      theme(legend.position="none", text = element_text(size=13)) +
      ggtitle(paste0("N = ", N, " , Simulation Setting = ", ID))+
      scale_fill_brewer(palette="BuPu")
  }
} 

fig1 = do.call(grid.arrange, p)
ggsave("results/figures/ME_compare.png", plot = fig1, width =  30, height = 35, units = "cm")



# PSR, NSR plots
p2<- list(); i = 0
for(N in c(100,400, 800, 1200)){
  for (ID in c(1,2,3)){
    i = i+1
    df  = read_delim(paste0("results/",N ,"_df", ID, "_gen.csv"), escape_double  =FALSE, delim = " ")
    df = as.tibble(df)
    p2[[i]] =  df %>% select('method', 'PSR', 'NSR') %>% 
      filter(method %in% c("ESL-LASSO(ours)", "PENSE-LASSO") )%>% 
      gather(key = "Type", value = "SR", -method) %>% 
      ggplot(aes(x = method , y = SR, fill = Type ))+
      geom_boxplot(width = .5)+ 
      theme( text = element_text(size=14) )+
      ggtitle(paste0("N = ", N, " , Simulation Setting = ", ID))+
      scale_fill_manual(values = wes_palette("Royal1", n = 2))
  }
} 

fig2 = do.call(grid.arrange, p2)
ggsave("results/figures/SR_compare.png", plot = fig2, width =  40, height = 35, units = "cm")







