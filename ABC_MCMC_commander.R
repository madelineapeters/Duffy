args = commandArgs(trailingOnly = TRUE)
x = as.integer(args[1])

mf.list = c(1/28470,0.001,0.01,0.1)
mf = mf.list[x]

## Save observed summary stats for LWK sample
obsSumStat = as.data.frame(matrix(nrow=1,ncol=16))
names(obsSumStat) = c("pi","SegSites","TajD","H","iHH","Fixed","Singletons","Doubletons","H1","H2","H12","H2.H1","haploNum","Frequency","EHH","Singletons/Fixed")
obsSumStat[1,] = c(1.17,14.00,-1.39,-5.52,3217.87,4985.00,5.00,2.00,0.57,0.02,0.74,0.03,12.00,1.00,0.67,0.0010)
write.csv(obsSumStat,"McManus_ObservedSumStat.csv",row.names=FALSE)

## AM1. Perform n simulations with parameters theta' randomly drawn from their priors, and each time compute associated set of summary statistics S'
source(paste(getwd(),"AM1_SLiM.r",sep="/"))

## AM2. Compute PLS components from the n theta' and S' vectors after a Box-Cox transformation of statistics
source(paste(getwd(),"AM2_SLiM.r",sep="/"))

## AM3. For all n simulations, transform the summary statistics S' into k retained PLS components, as S'_PLS. Transform the observed summary statistics S as S_PLS and compute P(mf|S_PLS,theta).
plot.opt = "yes" #Do you want to write a pdf file with the PLS-DA plot?
source(paste(getwd(),"AM3_SLiM.r",sep="/"))

## AM4. For all simulations for each model, calculate distance between observed and simulated summary statistics and retain the top 1% for each model. Then randomly sample 100 simulations (from the top 1% filtered) and for each calculate the mean distance from the drawn simulation to all the others (Euclidean distance of the PLS-DA transformed summary statistics). For each model, use Kolmogorov-Smirnov test to determine probability observed distance came from null distribution (distribution of pairwise Euclidean distances between drawn and all simulations for given model)
plot.opt = "yes" #Do you want to write a pdf file with the histogram plot?
source(paste(getwd(),"AM4_SLiM.r",sep="/"))
 