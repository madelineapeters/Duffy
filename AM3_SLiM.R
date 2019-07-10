AM1.BC.scaled = as.data.frame(AM1.pls$X)
comp.df = as.data.frame(AM1.pls$loadings$X)
comp.df$var = rownames(comp.df)

AM1.BC.trans = as.data.frame(t(sapply(1:nrow(AM1.BC.scaled),FUN=function(x){
  
  trans1 = sum(AM1.BC.scaled[x,]*comp.df$comp1)
  trans2 = sum(AM1.BC.scaled[x,]*comp.df$comp2)
  trans3 = sum(AM1.BC.scaled[x,]*comp.df$comp3)
  trans4 = sum(AM1.BC.scaled[x,]*comp.df$comp4)
  trans5 = sum(AM1.BC.scaled[x,]*comp.df$comp5)
  
  return(c(trans1,trans2,trans3,trans4,trans5))
  
})))
names(AM1.BC.trans) = paste("comp",1:5,sep="")
AM1.BC.trans$var = AM1.pls$Y
AM1.BC.trans$ID = 1:nrow(AM1.BC.trans)

# Centre and scale observed data, then transform in 5 PLS components
obsSumStat = read.csv(paste(getwd(),"McManus_ObservedSumStat.csv",sep="/")) %>% select(.,-iHH,-EHH,-pi)

obs.BC = t(sapply(3:ncol(AM1.df),FUN=function(x){
  if (any(AM1.df[,x]<=0)){
    yjPower(obsSumStat[1,x-2], lambda=p.list[x-2], jacobian.adjusted=FALSE)
  } else {
    bcPower(obsSumStat[1,x-2], lambda=p.list[x-2], jacobian.adjusted=FALSE)
  }
})) %>% as.data.frame()
names(obs.BC) = names(AM1.df)[3:ncol(AM1.df)]

center.constants = colMeans(AM1.BC[,1:(ncol(AM1.BC)-1)])
scale.sd = sapply(1:(ncol(AM1.BC)-1),FUN=function(x){
  sd(AM1.BC[,x])
})

obs.BC.scaled = obs.BC[1,] - center.constants
obs.BC.scaled[1,] = obs.BC.scaled[1,]/scale.sd

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
  
  AM1.BC.trans$Model = factor(AM1.BC.trans$var, levels=c("one","two","three","four"), labels=c("De Novo","0.1%","1%","10%"))
  
  ggplot()+geom_point(data=AM1.BC.trans,aes(x=comp1,y=comp2,col=Model))+geom_label(data=obs.trans,aes(x=comp1,y=comp2,label=var),size=6)+
    theme_classic()+theme(
      axis.title = element_text(size=15)
    )+
    xlab(paste("PLS component 1 (",100*round(AM1.pls$explained_variance$X[1],2),"%)",sep=""))+
    ylab(paste("PLS component 2 (",100*round(AM1.pls$explained_variance$X[2],2),"%)",sep=""))+
    scale_color_viridis(discrete=TRUE)
  
  ggsave("PLSDA_plot.tiff")
}
