# Load required packages, except mixOmics since it seems to upset dplyr
library(car)
library(pls)
library(dplyr)
library(EnvStats)

# Read in n simulations

AM1.1.df = read.csv(paste(getwd(),"AM1.1.csv",sep="/")) %>% 
  select(.,-seed,-EHH,-iHH,-pi) %>% 
  na.omit(.) %>% 
  mutate(.,mf = "one")

AM1.2.df = read.csv(paste(getwd(),"AM1.2.csv",sep="/")) %>% 
  select(.,-seed,-EHH,-iHH,-pi) %>% 
  na.omit(.) %>% 
  mutate(.,mf = "two")

AM1.3.df = read.csv(paste(getwd(),"AM1.3.csv",sep="/")) %>% 
  select(.,-seed,-EHH,-iHH,-pi) %>% 
  na.omit(.) %>% 
  mutate(.,mf = "three")

AM1.4.df = read.csv(paste(getwd(),"AM1.4.csv",sep="/")) %>% 
  select(.,-seed,-EHH,-iHH,-pi) %>% 
  na.omit(.) %>% 
  mutate(.,mf = "four")

AM1.df = bind_rows(AM1.1.df,AM1.2.df) %>% 
  bind_rows(.,AM1.3.df) %>% 
  bind_rows(.,AM1.4.df)

# Box-Cox transform summary statistics
BoxCoxTrans = function(y){
  lambda = 1
  constant = 0
  if(shapiro.test(AM1.df[,y])$p.value < 0.05) {
    if (min(AM1.df[,y]) > 0){
      lambda = boxcox(lm(AM1.df[,y]~AM1.df$s),optimize=TRUE)$lambda
    } else {
      lambda = powerTransform(AM1.df[,y]~AM1.df$s,family="yjPower")$lambda
    }
  }
  return(as.numeric(lambda))
}

p.list = sapply(3:ncol(AM1.df),FUN=BoxCoxTrans)

AM1.BC = as.data.frame(sapply(3:ncol(AM1.df),FUN=function(x){
  if (any(AM1.df[,x]<=0)){
    yjPower(AM1.df[,x], lambda=p.list[x-2], jacobian.adjusted=FALSE)
  } else {
    bcPower(AM1.df[,x], lambda=p.list[x-2], jacobian.adjusted=FALSE)
  }
}))
names(AM1.BC) = names(AM1.df)[3:ncol(AM1.df)]
AM1.BC = AM1.BC %>% bind_cols(.,select(AM1.df,mf))

# Compute PLS components from the n theta' and S' vectors (after Box-Cox transformation)
library(mixOmics)
AM1.pls = plsda(AM1.BC[,1:(ncol(AM1.BC)-1)],AM1.BC$mf,ncomp=5,scale=TRUE)
#plotIndiv(AM1.pls)
