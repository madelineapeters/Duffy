library(dplyr)
library(Zelig)

AM9.BC.trans = AM8.BC.trans

# Calculate distance between observed and simulated summary statistics for each of n simulations
Euclidean.fun = function(r){
  
  sqrt(sum(sapply(1:numPC,FUN=function(x){
    obs = obs.trans[1,x]
    sim = AM9.BC.trans[r,x]
    
    (obs - sim)^2
  })))
  
}
AM9.BC.trans$Distance = sapply(1:nrow(AM9.BC.trans),FUN=Euclidean.fun)

AM9.top = AM9.BC.trans %>% 
  top_frac(.,-0.1,Distance) %>% 
  mutate(.,choice = 1)

AM9.top$log.s = log10(AM9.top$s)

z.out = zelig(s ~ comp1 + comp2 + comp3 + comp4 + comp5, model = "normal.bayes", data=AM9.top, verbose=FALSE)
summary(z.out)

x.out = setx(z.out)
x.out = setx(z.out,comp1=obs.trans[1,1],comp2=obs.trans[1,2],comp3=obs.trans[1,3],comp4=obs.trans[1,4],comp5=obs.trans[1,5])
s.out1 = sim(z.out, x = x.out)
plot(s.out1)
s.df = zelig_qi_to_df(s.out1)

s.mean = mean(s.df$expected_value)
s.sd = sd(s.df$expected_value)

s.vals = as.data.frame(t(c(s.mean,s.sd)))
names(s.vals) = c("Mean","SD")

write.csv(s.vals,"sEstimate.csv",row.names=FALSE)

if (plot.opt == "yes"){
  
  prior.values = as.data.frame(10^runif(10000,-3,-0.5))
  names(prior.values) = "val"
  prior.values$Dist = "Prior"
  posterior.values = as.data.frame(s.df$expected_value)
  names(posterior.values) = "val"
  posterior.values$Dist = "Posterior"
  
  dist.df = bind_rows(prior.values,posterior.values)
  
  library(scales)
  
  ggplot(dist.df, aes(x=val, fill=Dist)) +
    geom_density(alpha=0.4) +
    theme_classic()+
    scale_fill_viridis(discrete=TRUE,option="E")+
    xlab("Selection coefficient")+ylab("Density")+
    theme_classic()+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    
    theme(
      axis.title = element_text(size=16),
      axis.text = element_text(size=12),
      legend.title = element_blank(),
      legend.text = element_text(size=12),
      legend.position = c(0.25,0.8)
    )
  p
  
  
  ggsave("PriorPost_dens.tiff")
  
}
  
  

