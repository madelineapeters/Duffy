source(paste(getwd(),"McManusCheck_SummaryStats.r",sep="/"))

runSLiM = function(x) {
  seed = x[1]
  s = x[2]
  mf = x[3]
  #cat("Running SLiM with seed ", seed, ", mu = ", mu, "\n");
  #output = system2("/plas1/apps/SLiM/bin/slim", c("-d", paste0("s=", s),"-d",paste0("mf=", mf),"-s", seed, " ~/Desktop/SLiM/McManusCheck.slim"), stdout=T)
  output = system2("/usr/local/bin/slim", c("-d", paste0("s=", s),"-d",paste0("mf=", mf),"-d",paste0("x=", mf),"-s", seed, " ~/SLiM/McManusCheck.slim"), stdout=T)
  #as.numeric(output[length(output)])
  
} 

s.prior = c(-3,-0.5)

n.sim = 10000

sum.stats = c("pi","SegSites","TajD","H","iHH","Fixed","Singletons","Doubletons","H1","H2","H12","H2.H1","haploNum","Frequency","EHH","Singletons/Fixed")

n.sim.df = as.data.frame(matrix(nrow=n.sim,ncol=length(sum.stats)+3))
names(n.sim.df) = c("mf","s",sum.stats,"seed")

c = 1
while (c < n.sim){
  
  print(paste("Working on simulation",c,sep=" "))
  
  tmp = try({
    seed = sample(1:100000, 1)
    exponent = runif(1,s.prior[1],s.prior[2])
    s = 10^exponent
    runSLiM(c(seed,s,mf))
    
    sum.stat.list = sum.stat.raw()  
    n.sim.df[c,] = c(mf,s,sum.stat.list,seed)
    
    write.csv(n.sim.df,paste("AM1",100*mf,"csv",sep="."),row.names=FALSE)
    
    c = c + 1
  })
  if(inherits(try, "try-error")) {
    recoverFromFailure()
  }
  
}
