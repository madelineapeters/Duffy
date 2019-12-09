library(dplyr)
library(Zelig)

AM5.top = AM3.BC.trans %>% 
  top_frac(.,-0.1,Distance) %>% 
  mutate(.,choice = 1)

AM5.top$Model = factor(AM5.top$var, levels=c("one","two","three","four","five"), labels=c(1,2,3,4,5))
AM5.top$Model = as.numeric(AM5.top$Model)

z.out = zelig(Model ~ comp1 + comp2 + comp3 + comp4 + comp5, model = "mlogit.bayes", data=AM5.top, verbose = FALSE)
summary(z.out)

x.out = setx(z.out,comp1=obs.trans[1,1],comp2=obs.trans[1,2],comp3=obs.trans[1,3],comp4=obs.trans[1,4],comp5=obs.trans[1,5])

s.out1 = sim(z.out, x = x.out)

s.df = zelig_qi_to_df(s.out1)

postProb = colMeans(s.df[7:11])
names(postProb) = c("De Novo","0.1%","1%","10%","25%")
BayesFactor = as.data.frame(outer(postProb, postProb, "/"))

write.csv(BayesFactor,"BayesFactor.csv",row.names=FALSE)
write.csv(postProb,"postProb.csv",row.names=FALSE)

