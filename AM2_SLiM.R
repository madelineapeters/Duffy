# Load required packages, except mixOmics since it seems to upset dplyr
library(car)
library(pls)
library(dplyr)
library(EnvStats)

# Read in n simulations

AM2.1.df = read.csv(paste(getwd(),"AM1.1.csv",sep="/")) %>% 
  slice(.,1:sims) %>% 
  select(.,-seed,-EHH,-iHH) %>% 
  na.omit(.) %>% 
  mutate(.,mf = "one")

AM2.2.df = read.csv(paste(getwd(),"AM1.2.csv",sep="/")) %>% 
  slice(.,1:sims) %>% 
  select(.,-seed,-EHH,-iHH) %>% 
  na.omit(.) %>% 
  mutate(.,mf = "two")

AM2.3.df = read.csv(paste(getwd(),"AM1.3.csv",sep="/")) %>% 
  slice(.,1:sims) %>% 
  select(.,-seed,-EHH,-iHH) %>% 
  na.omit(.) %>% 
  mutate(.,mf = "three")

AM2.4.df = read.csv(paste(getwd(),"AM1.4.csv",sep="/")) %>% 
  slice(.,1:sims) %>% 
  select(.,-seed,-EHH,-iHH) %>% 
  na.omit(.) %>% 
  mutate(.,mf = "four")

AM2.5.df = read.csv(paste(getwd(),"AM1.5.csv",sep="/")) %>% 
  slice(.,1:sims) %>% 
  select(.,-seed,-EHH,-iHH) %>% 
  na.omit(.) %>% 
  mutate(.,mf = "five")

AM2.df = bind_rows(AM2.1.df,AM2.2.df) %>% 
  bind_rows(.,AM2.3.df) %>% 
  bind_rows(.,AM2.4.df) %>% 
  bind_rows(.,AM2.5.df)

# Box-Cox transform summary statistics
BoxCoxTrans = function(y){
  lambda = 1
  constant = 0
  AM2.samp = sample_n(AM2.df,min(nrow(AM2.df),5000),replace=FALSE)
  if(shapiro.test(AM2.samp[,y])$p.value < 0.05) {
    if (min(AM2.df[,y]) > 0){
      lambda = EnvStats::boxcox(lm(AM2.df[,y]~AM2.df$s),optimize=TRUE)$lambda
    } else {
      lambda = powerTransform(AM2.df[,y]~AM2.df$s,family="yjPower")$lambda
    }
  }
  return(as.numeric(lambda))
}

p.list = sapply(3:ncol(AM2.df),FUN=BoxCoxTrans)

AM2.BC = as.data.frame(sapply(3:ncol(AM2.df),FUN=function(x){
  if (any(AM2.df[,x]<=0)){
    yjPower(AM2.df[,x], lambda=p.list[x-2], jacobian.adjusted=FALSE)
  } else {
    bcPower(AM2.df[,x], lambda=p.list[x-2], jacobian.adjusted=FALSE)
  }
}))
names(AM2.BC) = names(AM2.df)[3:ncol(AM2.df)]
AM2.BC = AM2.BC %>% bind_cols(.,select(AM2.df,mf))

# Compute PLS components from the n theta' and S' vectors (after Box-Cox transformation)
library(mixOmics)

AM2.pls = plsda(AM2.BC[,1:(ncol(AM2.BC)-1)],AM2.BC$mf,ncomp=5,scale=TRUE)
#plotIndiv(AM2.pls)
