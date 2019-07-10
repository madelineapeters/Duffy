# Calculate distance between observed and simulated summary statistics for each of n simulations
# Retain top 1% for each model (at 10% until more simulations done)
Euclidean.fun = function(r){
  
  sqrt(sum(sapply(1:5,FUN=function(x){
    obs = obs.trans[1,x]
    sim = AM3.BC.trans[r,x]
    
    (obs - sim)^2
  })))
  
}
AM3.BC.trans$Distance = sapply(1:nrow(AM3.BC.trans),FUN=Euclidean.fun)

AM4.top = AM3.BC.trans %>% group_by(var) %>% top_frac(.,-0.1,wt=Distance)

# For each model, randomly sample 100 simulations (from the top % filtered above) and for each calculate the mean distance from the drawn simulation to all the others (Euclidean distance of the PLS-DA transformed summary statistics)

AM4.top.1 = filter(AM4.top,var=="one")
r.list = sample(1:nrow(AM4.top.1),100,replace=FALSE)
distList.1 = sapply(r.list,FUN=function(r){
  
  AM4.comp = slice(AM4.top.1,-r)
  
  avgDist = mean(sapply(1:nrow(AM4.comp),FUN=function(x){
    
    comp = AM4.comp[x,1:5]
    sim = AM4.top.1[r,1:5]
    
    dist = sqrt(sum(as.numeric(c((comp - sim)^2))))
    
    return(dist)
    
  }))
  
  return(avgDist)
  
})

AM4.top.2 = filter(AM4.top,var=="two")
r.list = sample(1:nrow(AM4.top.2),100,replace=FALSE)
distList.2 = sapply(r.list,FUN=function(r){
  
  AM4.comp = slice(AM4.top.2,-r)
  
  avgDist = mean(sapply(1:nrow(AM4.comp),FUN=function(x){
    
    comp = AM4.comp[x,1:5]
    sim = AM4.top.2[r,1:5]
    
    dist = sqrt(sum(as.numeric(c((comp - sim)^2))))
    
    return(dist)
    
  }))
  
  return(avgDist)
  
})

AM4.top.3 = filter(AM4.top,var=="three")
r.list = sample(1:nrow(AM4.top.3),100,replace=FALSE)
distList.3 = sapply(r.list,FUN=function(r){
  
  AM4.comp = slice(AM4.top.3,-r)
  
  avgDist = mean(sapply(1:nrow(AM4.comp),FUN=function(x){
    
    comp = AM4.comp[x,1:5]
    sim = AM4.top.3[r,1:5]
    
    dist = sqrt(sum(as.numeric(c((comp - sim)^2))))
    
    return(dist)
    
  }))
  
  return(avgDist)
  
})

AM4.top.4 = filter(AM4.top,var=="four")
r.list = sample(1:nrow(AM4.top.4),100,replace=FALSE)
distList.4 = sapply(r.list,FUN=function(r){
  
  AM4.comp = slice(AM4.top.4,-r)
  
  avgDist = mean(sapply(1:nrow(AM4.comp),FUN=function(x){
    
    comp = AM4.comp[x,1:5]
    sim = AM4.top.4[r,1:5]
    
    dist = sqrt(sum(as.numeric(c((comp - sim)^2))))
    
    return(dist)
    
  }))
  
  return(avgDist)
  
})

# For each model, calculate the distance from the observed transformed summary statistics to all the simulated sets and calculated fraction of all simulated sets whose distance is below the observed distance

obsDist.1 = mean(AM4.top.1$Distance)
obsDist.2 = mean(AM4.top.2$Distance)
obsDist.3 = mean(AM4.top.3$Distance)
obsDist.4 = mean(AM4.top.4$Distance)

distList.1 = as.data.frame(distList.1)
names(distList.1) = "Distance"
distList.1$var = "De Novo"

distList.2 = as.data.frame(distList.2)
names(distList.2) = "Distance"
distList.2$var = "0.1%"

distList.3 = as.data.frame(distList.3)
names(distList.3) = "Distance"
distList.3$var = "1%"

distList.4 = as.data.frame(distList.4)
names(distList.4) = "Distance"
distList.4$var = "10%"

distList.df = bind_rows(distList.1,distList.2) %>% 
  bind_rows(.,distList.3) %>% 
  bind_rows(.,distList.4)

p1 = round(ks.test(obsDist.1,distList.1$Distance,alternative=c("two.sided"))$p.value,2)
p2 = round(ks.test(obsDist.2,distList.2$Distance,alternative=c("two.sided"))$p.value,2)
p3 = round(ks.test(obsDist.3,distList.3$Distance,alternative=c("two.sided"))$p.value,2)
p4 = round(ks.test(obsDist.4,distList.4$Distance,alternative=c("two.sided"))$p.value,2)

# Make dataframes for plotting histograms (in style of Figure 2 in McManus et al. Appendix)
obsDist.df = as.data.frame(matrix(nrow=4,ncol=3))
obsDist.df[,1] = c(obsDist.1,obsDist.2,obsDist.3,obsDist.4)
obsDist.df[,2] = c("De Novo","0.1%","1%","10%")
obsDist.df[,3] = c(p1,p2,p3,p4)
names(obsDist.df) = c("Distance","var","p")

# Plot histograms
if (plot.opt == "yes"){

  ggplot()+geom_histogram(data=distList.df,aes(Distance, stat(density)),fill=NA,col="black")+
    geom_vline(data=obsDist.df,aes(xintercept=Distance),col="dodgerblue",size=1)+
    geom_text(data=obsDist.df,aes(x=3.5,y=2,label=paste("p =",p,sep=" ")),size=4)+
    ylab("Density")+
    facet_wrap(var~.)+
    theme_classic()+
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size=14)
    )
  
  ggsave("KStest_hist.tiff")
  
  
}
