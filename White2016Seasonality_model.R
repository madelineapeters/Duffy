library(PBSddesolve)
library(tidyverse)

################################
# Functions
################################
SIV.fun = function (t, y,parms){
  
  #Set parameters
  ## Humans
  b = 0.5 #transmission probability (mosquito to human)
  r = 1/60 #rate of recovery from blood-stage infection
  delta = parms[1] #rate of recovery from dormancy
  gammaL = 1/223 #rate of latent hypnozoite clearance
  gammaD = 1/434 #rate of dormant hypnozoite clearance
  f = parms[2] #rate of relapse to blood stage infection
  if(t < 0){mu = 0} else {mu = parms[3]} #mortality rate due to blood stage infection
  h = parms[4] #dominance
  
  ## Mosquitoes
  a = parms[7] #mosquito biting frequency
  g = 0.1 #mosquito death rate
  n = 12 #duration of sporogony in mosquito
  c = 0.23 #transmission probability (human to mosquito)
  
  epsilon = 0.001 #low season density relative to mean
  kappa = parms[5] #seasonality parameter
  m0 = parms[6] #mean mosquito density over course of year
  t0 = 0.5 #time of peak density
  
  m = function(d){
    out = m0*(epsilon + (1-epsilon)*(base::pi/beta(0.5,0.5+kappa))*((1+cos(2*base::pi*((d/365)-0.5)))/2)^kappa)
    return(out)
  }
  
  #Variables to track (this is where R sets these equal to the initial values or updates them)
  
  #Duffy positive individuals
  S0p = y[1] #no blood stage infection, no latent stages, no dormant stages
  SLp = y[2] #no blood stage infection, has latent stages
  SDp = y[3] #no blood stage infection, has dormant stages
  SLDp = y[4] #no blood stage infection, has latent stages AND dormants stages
  
  I0p = y[5] #has blood stage infection, no latent stages, no dormant stages
  ILp = y[6] #has blood stage infection, has latent stages
  IDp = y[7] #has blood stage infection, has dormant stages
  ILDp = y[8] #has blood stage infection, has latent stages AND dormants stages
  
  #Duffy heterozygous individuals
  S0h = y[9] #no blood stage infection, no latent stages, no dormant stages
  SLh = y[10] #no blood stage infection, has latent stages
  SDh = y[11] #no blood stage infection, has dormant stages
  SLDh = y[12] #no blood stage infection, has latent stages AND dormants stages
  
  I0h = y[13] #has blood stage infection, no latent stages, no dormant stages
  ILh = y[14] #has blood stage infection, has latent stages
  IDh = y[15] #has blood stage infection, has dormant stages
  ILDh = y[16] #has blood stage infection, has latent stages AND dormants stages
  
  #Duffy negative individuals
  S0n = y[17]
  
  #Mosquitoes
  SM = y[18]
  IM = y[19]
  
  #force of infection in humans
  lambda = m(t)*a*b*IM
  
  #Dynamics of Duffy positive individuals
  dS0pdt = -lambda*S0p + r*I0p + gammaL*SLp + gammaD*SDp #done
  dSLpdt = - (lambda + f + gammaL)*SLp + r*ILp + (delta + gammaD)*SLDp + delta*SDp #done
  dSDpdt = - (lambda + gammaD + delta)*SDp + r*IDp + gammaL*SLDp #done
  dSLDpdt = - (lambda + f + gammaD + delta + gammaL)*SLDp + r*ILDp #done
  
  dI0pdt = - (lambda + r)*I0p + gammaL*ILp + gammaD*IDp - mu*I0p #done
  dILpdt = f*SLp - (r + lambda + gammaL)*ILp + (gammaD + delta)*ILDp + delta*IDp - mu*ILp #done
  dIDpdt = lambda*(SDp + S0p + I0p) - (r + gammaD + delta)*IDp + gammaL*ILDp - mu*IDp #done
  dILDpdt = lambda*(SLp + ILp + SLDp) + f*SLDp - (gammaD + delta + r + gammaL)*ILDp - mu*ILDp #done
  
  #Dynamics of Duffy heterozygous individuals
  dS0hdt = -lambda*h*S0h + r*I0h + gammaL*SLh + gammaD*SDh #done
  dSLhdt = - (lambda*h + f + gammaL)*SLh + r*ILh + (delta + gammaD)*SLDh + delta*SDh #done
  dSDhdt = - (lambda*h + gammaD + delta)*SDh + r*IDh + gammaL*SLDh #done
  dSLDhdt = - (lambda*h + f + gammaD + delta + gammaL)*SLDh + r*ILDh #done
  
  dI0hdt = - (lambda*h + r)*I0h + gammaL*ILh + gammaD*IDh - mu*I0h
  dILhdt = f*SLh - (r + lambda*h + gammaL)*ILh + (gammaD + delta)*ILDh + delta*IDh - mu*ILh
  dIDhdt = lambda*h*(SDh + S0h + I0h) - (r + gammaD + delta)*IDh + gammaL*ILDh - mu*IDh
  dILDhdt = lambda*h*(SLh + ILh + SLDh) + f*SLDh - (gammaD + delta + r + gammaL)*ILDh - mu*ILDh
  
  #Dynamics of Duffy negative individuals
  dS0ndt = 0
  
  #Dynamics of mosquitoes
  dSMdt = g - a*c*(I0p + ILp + IDp + ILDp + I0h + ILh + IDh +ILDh)*(exp(-g*n) - IM) - g*SM 
  dIMdt = a*c*(I0p + ILp + IDp + ILDp + I0h + ILh + IDh +ILDh)*(exp(-g*n) - IM) - g*IM
  
  return(c(dS0pdt,dSLpdt,dSDpdt,dSLDpdt,dI0pdt,dILpdt,dIDpdt,dILDpdt,dS0hdt,dSLhdt,dSDhdt,dSLDhdt,dI0hdt,dILhdt,dIDhdt,dILDhdt,dS0ndt,dSMdt,dIMdt))
  
}
SIV.fun.nodorm = function (t, y,parms){
  #Set parameters
  ## Humans
  b = 0.5 #transmission probability (mosquito to human)
  r = 1/60 #rate of recovery from blood-stage infection
  gammaL = 1/223 #rate of latent hypnozoite clearance
  f = parms[1] #rate of relapse to blood stage infection
  if(t < 0){mu = 0} else {mu = parms[2]} #mortality rate due to blood stage infection
  h = parms[3] #dominance
  
  ## Mosquitoes
  a = parms[6] #mosquito biting frequency
  g = 0.1 #mosquito death rate
  n = 12 #duration of sporogony in mosquito
  c = 0.23 #transmission probability (human to mosquito)
  
  epsilon = 0.001 #low season density relative to mean
  kappa = parms[4] #seasonality parameter
  m0 = parms[5] #mean mosquito density over course of year
  t0 = 0.5 #time of peak density
  
  m = function(d){
    out = m0*(epsilon + (1-epsilon)*(base::pi/beta(0.5,0.5+kappa))*((1+cos(2*base::pi*((d/365)-0.5)))/2)^kappa)
    return(out)
  }
  
  #Variables to track (this is where R sets these equal to the initial values or updates them)
  
  #Duffy positive individuals
  S0p = y[1] #no blood stage infection, no latent stages, no dormant stages
  SLp = y[2] #no blood stage infection, has latent stages
  
  I0p = y[3] #has blood stage infection, no latent stages, no dormant stages
  ILp = y[4] #has blood stage infection, has latent stages
  
  #Duffy heterozygous individuals
  S0h = y[5] #no blood stage infection, no latent stages, no dormant stages
  SLh = y[6] #no blood stage infection, has latent stages
  
  I0h = y[7] #has blood stage infection, no latent stages, no dormant stages
  ILh = y[8] #has blood stage infection, has latent stages
  
  #Duffy negative individuals
  S0n = y[9]
  
  #Mosquitoes
  SM = y[10]
  IM = y[11]
  
  #force of infection in humans
  lambda = m(t)*a*b*IM
  
  #Dynamics of Duffy positive individuals
  dS0pdt = -lambda*S0p + r*I0p + gammaL*SLp #done
  dSLpdt = - (lambda + f + gammaL)*SLp + r*ILp #done
  
  dI0pdt = - (lambda + r)*I0p + gammaL*ILp - mu*I0p #done
  dILpdt = lambda*(S0p + SLp + I0p) + f*SLp - (r + gammaL)*ILp - mu*ILp #done
  
  #Dynamics of Duffy heterozygous individuals
  dS0hdt = -lambda*h*S0h + r*I0h + gammaL*SLh #done
  dSLhdt = - (lambda*h + f + gammaL)*SLh + r*ILh #done
  
  dI0hdt = - (lambda*h + r)*I0h + gammaL*ILh - mu*I0h
  dILhdt = lambda*h*(S0h + SLh + I0h) + f*SLh - (r + gammaL)*ILh - mu*ILh #done
  
  #Dynamics of Duffy negative individuals
  dS0ndt = 0
  
  #Dynamics of mosquitoes
  dSMdt = g - a*c*(I0p + ILp + I0h + ILh)*(exp(-g*n) - IM) - g*SM 
  dIMdt = a*c*(I0p + ILp + I0h + ILh)*(exp(-g*n) - IM) - g*IM
  
  return(c(dS0pdt,dSLpdt,dI0pdt,dILpdt,dS0hdt,dSLhdt,dI0hdt,dILhdt,dS0ndt,dSMdt,dIMdt))
  
}

################################
# Set-up
################################
#Model parameters
# delta = rate of recovery from dormancy (-> infinity for no dormancy) (1/162)
delta = 1/162
# f = rate of relapse (-> 0 for no relapse) (1/72)
f = 1/72
# mu = mortality due to blood stage infection
## Mean duration of infection is 60 days
## For mortality rate of 0.00012 (0.12/1000 infections)
mu = 1/500000

#Generation time
gen.time = 29

#Population size at start of generation
pop.size = 1000

################################
# Run 
################################    
#Interested in changing mean density of mosquitoes m0, seasonality parameter kappa, dominance h
p0.list = c(0.99,seq(0.95,0.05,-0.05),0.04,0.03,0.02,0.01) #frequency of Duffy positive allele at start of generation
m.list = c(0.05,0.5,1,2,5) #average annual mosquito density
kappa.list = c(0,0.37,1) #seasonality of mosquitoes
h.list = c(0,0.5,1) #dominance of Duffy positive allele (h = 0.5 means heterozygote 50% less susceptible to infection than Duffy positive homozygote; h = 0 means heterozygotes complete protected)
a.list = c(0.05,0.1,0.21,0.5,1)

run.outputs = as.data.frame(matrix(nrow=1,ncol=10))
names(run.outputs) = c("p0","m0","kappa","h","a","Positive.survival","Heterozygous.survival","Negative.survival","Positive.selection","Heterozygous.selection")

iteration = 1

comment.out = function(){for (p0.index in 1:24){
  for (m.index in 1:5){
    for (kappa.index in 1:3){
      for (h.index in 1:3){
        for (a.index in 1:5){
          print(paste("Working on iteration",iteration,sep=" "))
          
          p.set = p0.list[p0.index]
          h = h.list[h.index]
          kappa = kappa.list[kappa.index]
          m0 = m.list[m.index]
          a = a.list[a.index]
          
          #Set initial numbers of Duffy positive, heterozygous and negative individuals
          Positive.start = pop.size*p.set*p.set
          Heterozygous.start = pop.size*2*p.set*(1-p.set)
          Negative.start = pop.size*(1-p.set)*(1-p.set)
          
          yinit = c(y1 = Positive.start*0.5, #S0p
                    y2 = Positive.start*0.5, #SLp
                    y3 = 0, #SDp
                    y4 = 0, #SLDp
                    y5 = 0, #I0p
                    y6 = 0, #ILp
                    y7 = 0, #IDp
                    y8 = 0, #ILDp
                    y9 = Heterozygous.start*(0.5+0.5*(1-h)), #S0h
                    y10 = Heterozygous.start*0.5*h, #SLh
                    y11 = 0, #SDh
                    y12 = 0, #SLDh
                    y13 = 0, #I0h
                    y14 = 0, #ILh
                    y15 = 0, #IDh
                    y16 = 0, #ILDH
                    y17 = Negative.start, #S0N
                    y18 = 0, #SM
                    y19 = 1) #IM
          times = -3000:(gen.time*365)
          
          parameters = c(delta,f,mu,h,kappa,m0,a)
          
          ###Run DDE
          dde.SIV = dde(y = yinit, times = times, func = SIV.fun, parms = parameters, hbsize = 50000)
          names(dde.SIV) = c("Time",
                             paste(c("S0","SL","SD","SLD","I0","IL","ID","ILD"),"p",sep="."),
                             paste(c("S0","SL","SD","SLD","I0","IL","ID","ILD"),"h",sep="."),
                             "S0n","SM","IM")
          
          dde.SIV = dde.SIV %>% 
            mutate(.,"Positive" = S0.p + SL.p + SLD.p + SD.p + I0.p + IL.p + ID.p + ILD.p) %>% 
            mutate(.,"Heterozygous" = S0.h + SL.h + SLD.h + SD.h + I0.h + IL.h + ID.h + ILD.h) %>% 
            mutate(.,"Negative" = S0n) %>% 
            mutate(.,"Mosquitoes" = SM + IM) %>% 
            mutate(.,"BloodStage" = I0.p + IL.p + ID.p + ILD.p + I0.h + IL.h + ID.h + ILD.h) %>% 
            mutate(.,"Latent" = SL.p + SLD.p + IL.p + ILD.p + SL.h + SLD.h + IL.h + ILD.h) %>% 
            mutate(.,"Dormant" = SD.p + SLD.p + ID.p + ILD.p + SD.h + SLD.h + ID.h + ILD.h) %>% 
            mutate(.,"p" = (Positive + 0.5*Heterozygous)/(Positive + Heterozygous + Negative))
          
          final.p = dde.SIV$p[nrow(dde.SIV)]
          final.Positive = dde.SIV$Positive[nrow(dde.SIV)]
          final.Heterozygous = dde.SIV$Heterozygous[nrow(dde.SIV)]
          final.Negative = dde.SIV$Negative[nrow(dde.SIV)]
          
          #Calculate survival rates
          Positive.survival = final.Positive/Positive.start
          Heterozygous.survival = final.Heterozygous/Heterozygous.start
          Negative.survival = final.Negative/Negative.start
          
          #Calculate fitness values by dividing by maximum survival rate (should be Duffy negative survival, which should be 1, if Duffy negative individuals are entirely protected and there is no cost)
          max.survival = max(c(Positive.survival,Heterozygous.survival,Negative.survival))
          Positive.fitness = Positive.survival/max.survival
          Heterozygous.fitness = Heterozygous.survival/max.survival
          Negative.fitness = Negative.survival/max.survival
          
          #Calculate selection coefficients, which are just 1 - fitness
          Positive.selection = 1 - Positive.survival
          Heterozygous.selection = 1 - Heterozygous.survival
          
          run.outputs[1,] = c(p.set,m0,kappa,h,a,Positive.fitness,Heterozygous.fitness,Negative.fitness,Positive.selection,Heterozygous.selection)
          if (p0.index == 1 & m.index == 1 & kappa.index == 1 & h.index == 1 & a.index == 1){
            joint.outputs = run.outputs
          } else {
            joint.outputs = bind_rows(joint.outputs,run.outputs)
          }
          
          iteration = iteration + 1
          
          write.csv(joint.outputs,"White2016Seasonality_outputs.csv",row.names=FALSE)
        }
      }
    }
  }
}}

iteration = 1

for (p0.index in 1:24){
  for (m.index in 1:5){
    for (kappa.index in 1:3){
      for (h.index in 1:3){
        for (a.index in 1:5){
          print(paste("Working on iteration",iteration,sep=" "))
          
          p.set = p0.list[p0.index]
          h = h.list[h.index]
          kappa = kappa.list[kappa.index]
          m0 = m.list[m.index]
          a = a.list[a.index]
          
          #Set initial numbers of Duffy positive, heterozygous and negative individuals
          Positive.start = pop.size*p.set*p.set
          Heterozygous.start = pop.size*2*p.set*(1-p.set)
          Negative.start = pop.size*(1-p.set)*(1-p.set)
          
          yinit = c(y1 = Positive.start*0.5, #S0p
                    y2 = Positive.start*0.5, #SLp
                    y3 = 0, #I0p
                    y4 = 0, #ILp
                    y5 = Heterozygous.start*(0.5+0.5*(1-h)), #S0h
                    y6 = Heterozygous.start*0.5*h, #SLh
                    y7 = 0, #I0h
                    y8 = 0, #ILh
                    y9 = Negative.start, #S0N
                    y10 = 0, #SM
                    y11 = 1) #IM
          times = -3000:(gen.time*365)
          
          parameters = c(f,mu,h,kappa,m0,a)
          
          ###Run DDE
          dde.SIV = dde(y = yinit, times = times, func = SIV.fun.nodorm, parms = parameters, hbsize = 50000)
          names(dde.SIV) = c("Time",
                             paste(c("S0","SL","I0","IL"),"p",sep="."),
                             paste(c("S0","SL","I0","IL"),"h",sep="."),
                             "S0n","SM","IM")
          
          dde.SIV = dde.SIV %>% 
            mutate(.,"Positive" = S0.p + SL.p + I0.p + IL.p) %>% 
            mutate(.,"Heterozygous" = S0.h + SL.h + I0.h + IL.h) %>% 
            mutate(.,"Negative" = S0n) %>% 
            mutate(.,"Mosquitoes" = SM + IM) %>% 
            mutate(.,"BloodStage" = I0.p + IL.p + I0.h + IL.h) %>% 
            mutate(.,"Latent" = SL.p + IL.p + SL.h + IL.h) %>% 
            mutate(.,"p" = (Positive + 0.5*Heterozygous)/(Positive + Heterozygous + Negative))
          
          final.p = dde.SIV$p[nrow(dde.SIV)]
          final.Positive = dde.SIV$Positive[nrow(dde.SIV)]
          final.Heterozygous = dde.SIV$Heterozygous[nrow(dde.SIV)]
          final.Negative = dde.SIV$Negative[nrow(dde.SIV)]
          
          #Calculate survival rates
          Positive.survival = final.Positive/Positive.start
          Heterozygous.survival = final.Heterozygous/Heterozygous.start
          Negative.survival = final.Negative/Negative.start
          
          #Calculate fitness values by dividing by maximum survival rate (should be Duffy negative survival, which should be 1, if Duffy negative individuals are entirely protected and there is no cost)
          max.survival = max(c(Positive.survival,Heterozygous.survival,Negative.survival))
          Positive.fitness = Positive.survival/max.survival
          Heterozygous.fitness = Heterozygous.survival/max.survival
          Negative.fitness = Negative.survival/max.survival
          
          #Calculate selection coefficients, which are just 1 - fitness
          Positive.selection = 1 - Positive.survival
          Heterozygous.selection = 1 - Heterozygous.survival
          
          run.outputs[1,] = c(p.set,m0,kappa,h,a,Positive.fitness,Heterozygous.fitness,Negative.fitness,Positive.selection,Heterozygous.selection)
          if (p0.index == 1 & m.index == 1 & kappa.index == 1 & h.index == 1 & a.index == 1){
            joint.outputs = run.outputs
          } else {
            joint.outputs = bind_rows(joint.outputs,run.outputs)
          }
          
          iteration = iteration + 1
          
          write.csv(joint.outputs,"White2016Seasonality_outputs_nodorm.csv",row.names=FALSE)
        }
      }
    }
  }
}



  