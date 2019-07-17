# Load required packages, except mixOmics since it seems to upset dplyr
library(car)
library(pls)
library(dplyr)
library(EnvStats)

# Read in n simulations

AM7.df = read.csv(paste(getwd(),"AM6.csv",sep="/")) %>% 
  slice(.,1:sims) %>% 
  select(.,-seed,-EHH,-iHH) %>% 
  na.omit(.)

# Box-Cox transform summary statistics
BoxCoxTrans = function(y){
  lambda = 1
  constant = 0
  AM7.samp = sample_n(AM7.df,min(nrow(AM7.df),5000),replace=FALSE)
  if(shapiro.test(AM7.samp[,y])$p.value < 0.05) {
    if (min(AM7.df[,y]) > 0){
      lambda = EnvStats::boxcox(lm(AM7.df[,y]~AM7.df$s),optimize=TRUE)$lambda
    } else {
      lambda = powerTransform(AM7.df[,y]~AM7.df$s,family="yjPower")$lambda
    }
  }
  return(as.numeric(lambda))
}

p.list = sapply(3:ncol(AM7.df),FUN=BoxCoxTrans)

AM7.BC = as.data.frame(sapply(3:ncol(AM7.df),FUN=function(x){
  if (any(AM7.df[,x]<=0)){
    yjPower(AM7.df[,x], lambda=p.list[x-2], jacobian.adjusted=FALSE)
  } else {
    bcPower(AM7.df[,x], lambda=p.list[x-2], jacobian.adjusted=FALSE)
  }
}))
names(AM7.BC) = names(AM7.df)[3:ncol(AM7.df)]
AM7.BC = AM7.BC %>% bind_cols(.,select(AM7.df,mf))

# Compute PCA components from the n theta' and S' vectors (after Box-Cox transformation)
AM7.pca = prcomp(AM7.BC[,1:(ncol(AM7.BC)-1)],scale=TRUE)
#plotIndiv(AM7.pls)

# Compute number of principle components to return to capture 95% of variance
numPC = min(which(cumsum((AM7.pca$sdev^2)/sum(AM7.pca$sdev^2)) > 0.95))
