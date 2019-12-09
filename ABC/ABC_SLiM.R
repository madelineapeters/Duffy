library(EasyABC)

runSLiM = function(x) {
  seed = x[1]
  s = x[2]
  #cat("Running SLiM with seed ", seed, ", mu = ", mu, "\n");
  output = system2("/usr/local/bin/slim", c("-d", paste0("s=", s),
                                             "-s", seed, " ~/Desktop/SLiM/ABC.slim"), stdout=T)
  as.numeric(output[length(output)])
  
} 

prior = list(c("unif",0.002,0.07))
observed = 0.98

#ABC_SLiM = ABC_sequential(method="Lenormand", use_seed=TRUE, model=runSLiM, prior=prior, summary_stat_target=observed, nb_simul=1000)

ABC_SLiM_MCMC = ABC_mcmc(method="Marjoram", use_seed=TRUE, model=runSLiM, prior=prior, summary_stat_target=observed, n_cluster = 10)


x = sum(ABC_SLiM$param * ABC_SLiM$weights)

print(x)
