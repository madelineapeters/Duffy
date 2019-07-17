detach(package:mixOmics)
detach(package:dplyr)
library(dplyr)

AM8.BC.trans = as.data.frame(AM7.pca$x) %>% select(.,1:numPC)
names(AM8.BC.trans) = paste("comp",1:numPC,sep="")
AM8.BC.trans$s = AM7.df$s
AM8.BC.trans$ID = 1:nrow(AM8.BC.trans)

# Centre and scale observed data, then transform in 5 PLS components
obsSumStat = read.csv(paste(getwd(),"McManus_ObservedSumStat.csv",sep="/")) %>% select(.,-iHH,-EHH)

obs.BC = t(sapply(3:ncol(AM7.df),FUN=function(x){
  if (any(AM7.df[,x]<=0)){
    yjPower(obsSumStat[1,x-2], lambda=p.list[x-2], jacobian.adjusted=FALSE)
  } else {
    bcPower(obsSumStat[1,x-2], lambda=p.list[x-2], jacobian.adjusted=FALSE)
  }
})) %>% as.data.frame()
names(obs.BC) = names(AM7.df)[3:ncol(AM7.df)]

center.constants = colMeans(AM7.BC[,1:(ncol(AM7.BC)-1)])
scale.sd = sapply(1:(ncol(AM7.BC)-1),FUN=function(x){
  sd(AM7.BC[,x])
})

obs.BC.scaled = obs.BC[1,] - center.constants
obs.BC.scaled[1,] = obs.BC.scaled[1,]/scale.sd

comp.df = as.data.frame(AM7.pca$rotation[,1:numPC])
names(comp.df) = paste("comp",1:numPC,sep="")

obs.trans = as.data.frame(t(sapply(1,FUN=function(x){
  
  trans1 = sum(obs.BC.scaled[x,]*comp.df$comp1)
  trans2 = sum(obs.BC.scaled[x,]*comp.df$comp2)
  trans3 = sum(obs.BC.scaled[x,]*comp.df$comp3)
  trans4 = sum(obs.BC.scaled[x,]*comp.df$comp4)
  trans5 = sum(obs.BC.scaled[x,]*comp.df$comp5)
  
  return(c(trans1,trans2,trans3,trans4,trans5))
  
})))
names(obs.trans) = paste("comp",1:5,sep="")
obs.trans$var = "Observed"

if (plot.opt == "yes"){
  library(ggplot2)
  library(viridis)
  
  ggplot()+geom_point(data=AM8.BC.trans,aes(x=comp1,y=comp2,col=s))+geom_label(data=obs.trans,aes(x=comp1,y=comp2,label=var),size=6)+
    theme_classic()+theme(
      axis.title = element_text(size=15)
    )+
    xlab(paste("PC 1 (",100*round(((AM7.pca$sdev^2)/sum(AM7.pca$sdev^2))[1],2),"%)",sep=""))+
    ylab(paste("PC 2 (",100*round(((AM7.pca$sdev^2)/sum(AM7.pca$sdev^2))[2],2),"%)",sep=""))+
    scale_color_viridis()
  
  ggsave("PCA_plot.tiff")
}
