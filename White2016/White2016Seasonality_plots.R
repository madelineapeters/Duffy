library(tidyverse)
library(ggplot2)
library(gridExtra)

joint.outputs = read.csv("~/Desktop/Mideo.lab2/Duffy/Transmission_models/White2016/White2016Seasonality_outputs.csv")

joint.outputs$m0 = factor(joint.outputs$m0,levels=c(0.05,0.1,0.5,1.0,2.0), labels=c("m[0] == 0.05","m[0] == 0.1", "m[0] == 0.5", "m[0] == 1.0","m[0] == 2.0"))

selection.facet1 = ggplot()+
  geom_line(data=filter(joint.outputs,m0%in%c("m[0] == 0.5","m[0] == 1.0","m[0] == 2.0"),h%in%c(0,0.5,1.0),a==0.05),aes(x=1-p0,y=Positive.selection,col=factor(h)))+
  geom_line(data=filter(joint.outputs,m0%in%c("m[0] == 0.5","m[0] == 1.0","m[0] == 2.0"),h%in%c(0,0.5,1.0),a==0.05),aes(x=1-p0,y=Heterozygous.selection,col=factor(h)),linetype="dashed")+
  facet_grid(kappa~m0,labeller=label_parsed)+
  ylim(-0.0001,0.02)+
  theme_bw()+
  scale_x_continuous(breaks=c(0,0.5,1.0),labels=c("0","0.5","1.0"))+
  #scale_y_continuous(breaks=c(0,0.0075,0.015),labels=c("0","7.5e-3","1.5e-2"))+
  scale_color_manual(name="Dominance",values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"))+
  xlab("Frequency of Duffy null allele")+
  ylab("Selection coefficient")+
  ggtitle("With dormancy (d = 162), a = 0.05")+
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    strip.text.y = element_blank(),
    legend.position="left"
  )

selection.facet2 = ggplot()+
  geom_line(data=filter(joint.outputs,m0%in%c("m[0] == 0.5","m[0] == 1.0","m[0] == 2.0"),h%in%c(0,0.5,1.0),a==0.21),aes(x=1-p0,y=Positive.selection,col=factor(h)))+
  geom_line(data=filter(joint.outputs,m0%in%c("m[0] == 0.5","m[0] == 1.0","m[0] == 2.0"),h%in%c(0,0.5,1.0),a==0.21),aes(x=1-p0,y=Heterozygous.selection,col=factor(h)),linetype="dashed")+
  facet_grid(kappa~m0,labeller=label_parsed)+
  theme_bw()+
  ylim(-0.0001,0.02)+
  scale_x_continuous(breaks=c(0,0.5,1.0),labels=c("0","0.5","1.0"))+
  #scale_y_continuous(breaks=c(0,0.0075,0.015),labels=c("0","7.5e-3","1.5e-2"))+
  scale_color_manual(name="Dominance",values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"))+
  ggtitle("a = 0.21")+
  xlab("Frequency of Duffy null allele")+
  ylab("Selection coefficient")+
  ggtitle("With dormancy (d = 162), a = 0.21")+
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    strip.text.y = element_blank(),
    legend.position="none"
  )

kappa.df = data.frame()
kappa.list = c(0,0.37,1)
for (kappa in kappa.list){
  m0 = 1
  m = function(d){
    out = m0*(epsilon + (1-epsilon)*(base::pi/beta(0.5,0.5+kappa))*((1+cos(2*base::pi*((d/365)-0.5)))/2)^kappa)
    return(out)
  }
  kappa.temp = as.data.frame(sapply(seq(0,365,1),m))
  names(kappa.temp) = "Density"
  kappa.temp$kappa = kappa
  kappa.temp$m0 = "Seasonality"
  kappa.temp$Time = seq(0,365,1)
  kappa.df = bind_rows(kappa.df,kappa.temp)
  
}

kappa.facet = ggplot(data=kappa.df)+geom_line(data=kappa.df,aes(x=Time,y=Density))+
  facet_grid(kappa~m0,labeller = label_bquote(kappa == .(kappa)))+
  theme_bw()+
  scale_x_continuous(breaks=cumsum(c(0,31,28,31,30,31,30,31,31,30,31,30,31)),labels=0:12)+
  xlab("Month")+
  ylab("Density")+
  ggtitle("")+
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank()#,
    #strip.text.x = element_text(colour = 'white')
  )


###################
### No dormancy ###
###################
joint.outputs = read.csv("~/Desktop/Mideo.lab2/Duffy/Transmission_models/White2016/White2016Seasonality_outputs_nodorm.csv")

joint.outputs$m0 = factor(joint.outputs$m0,levels=c(0.05,0.1,0.5,1.0,2.0), labels=c("m[0] == 0.05","m[0] == 0.1", "m[0] == 0.5", "m[0] == 1.0","m[0] == 2.0"))

selection.facet3 = ggplot()+
  geom_line(data=filter(joint.outputs,m0%in%c("m[0] == 0.5","m[0] == 1.0","m[0] == 2.0"),h%in%c(0,0.5,1.0),a==0.05),aes(x=1-p0,y=Positive.selection,col=factor(h)))+
  geom_line(data=filter(joint.outputs,m0%in%c("m[0] == 0.5","m[0] == 1.0","m[0] == 2.0"),h%in%c(0,0.5,1.0),a==0.05),aes(x=1-p0,y=Heterozygous.selection,col=factor(h)),linetype="dashed")+
  facet_grid(kappa~m0,labeller=label_parsed)+
  theme_bw()+
  ylim(-0.001,0.0225)+
  scale_x_continuous(breaks=c(0,0.5,1.0),labels=c("0","0.5","1.0"))+
  #scale_y_continuous(breaks=c(0,0.0075,0.015),labels=c("0","7.5e-3","1.5e-2"))+
  scale_color_manual(name="Dominance",values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"))+
  xlab("Frequency of Duffy null allele")+
  ylab("Selection coefficient")+
  ggtitle("No dormancy, a = 0.05")+
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    strip.text.y = element_blank(),
    legend.position="left"
  )
#selection.facet1

selection.facet4 = ggplot()+
  geom_line(data=filter(joint.outputs,m0%in%c("m[0] == 0.5","m[0] == 1.0","m[0] == 2.0"),h%in%c(0,0.5,1.0),a==0.21),aes(x=1-p0,y=Positive.selection,col=factor(h)))+
  geom_line(data=filter(joint.outputs,m0%in%c("m[0] == 0.5","m[0] == 1.0","m[0] == 2.0"),h%in%c(0,0.5,1.0),a==0.21),aes(x=1-p0,y=Heterozygous.selection,col=factor(h)),linetype="dashed")+
  facet_grid(kappa~m0,labeller=label_parsed)+
  theme_bw()+
  ylim(-0.001,0.0225)+
  scale_x_continuous(breaks=c(0,0.5,1.0),labels=c("0","0.5","1.0"))+
  #scale_y_continuous(breaks=c(0,0.0075,0.015),labels=c("0","7.5e-3","1.5e-2"))+
  scale_color_manual(name="Dominance",values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"))+
  ggtitle("No dormancy, a = 0.21")+
  xlab("Frequency of Duffy null allele")+
  ylab("Selection coefficient")+
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    strip.text.y = element_blank(),
    legend.position="none"
  )
#selection.facet2

grid.arrange(selection.facet1,selection.facet2,kappa.facet,
             selection.facet3,selection.facet4,kappa.facet,
             nrow=2,ncol=3,widths=c(2.4,2,1.25))



