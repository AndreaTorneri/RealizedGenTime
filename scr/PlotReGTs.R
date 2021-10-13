###########################################################
###########################################################
## Plot Fig 1 - Contact Rates
###########################################################
###########################################################
library(ggplot2)
library(latex2exp)
library(ggpubr)

f1<-function(t){
  value<-dweibull(t,shape = 3.5, scale = 0.56)/(5*0.9995042)
  return(value)
}

f2<-function(t){
  a<-exp(-0.5*(9*t))
  d<-2*(1-exp(-9))
  int.v<-5*0.1098903
  return((a/d)/int.v)
}

f3<-function(t){
  d<-5
  value<-1/d
  return(value)
}
t<-seq(0,1,0.001)

im1<-f1(t)*5*1.5
im2<-f2(t)*5*1.5
im3<-f3(t)*5*1.5

ds1<-data.frame(X=t, Rate1=im3, Rate2=im2, Rate3=im1)
ds1$IPD<-factor(dfR0$IPD, levels=c("Const(1)","Gamma(2,0.5)","Exp(1)","Gamma(0.5,2)"))




colors <- c("Peaked" = "#E41A1C", "Decreasing" = "#377EB8", "Constant" = "#4DAF4A")
ggplot(ds1, aes(x = X)) +
  geom_line(aes(y = im1, color = "Peaked"), size = 1) +
  geom_line(aes(y = im2, color = "Decreasing"), size = 1) +
  geom_line(aes(y = im3, color = "Constant"), size = 1) +
  labs(x = "Time", y="Effective Contact Profile Value",
       color = "Profiles:") +
  scale_color_manual(values = colors, guide = guide_legend(reverse = TRUE))+theme(
    #axis.title.y =element_blank(), 
    #axis.text.y = element_blank(), 
    #axis.ticks.y = element_blank(),
    
    panel.background = element_rect(fill = "white",                                                                                                                       
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                    colour = "gray"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "white"),  
    legend.background = element_rect(fill = "white"),  
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.size = unit(2, "cm"),
    legend.key.width = unit(2,"cm"),
    legend.title = element_text(size=12, face = "bold"),
    legend.text = element_text(size=11),
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=11),
    axis.title = element_text(size=12),
  )+ylim(c(0,6.5))

#plot(a3)

######################################################################
# VIOLIN PLOTS
######################################################################
library(ggplot2)
library(ggpubr)
library(Hmisc)

######################################################################
#### R0
######################################################################
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R1.5_mu1_IPDistConst_RateConst_.RData")
MeanGTR1.5<-MeanGT
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R2_mu1_IPDistConst_RateConst_.RData")
MeanGTR2<-MeanGT
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistConst_RateConst_.RData")
MeanGTR3<-MeanGT
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R4_mu1_IPDistConst_RateConst_.RData")
MeanGTR4<-MeanGT
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R5_mu1_IPDistConst_RateConst_.RData")
MeanGTR5<-MeanGT

dfR0<-data.frame(MeanRGT=c(MeanGTR1.5,MeanGTR2,MeanGTR3,MeanGTR4, MeanGTR5),R0= c(rep("R0=1.5",length(MeanGTR1.5)),rep("R0=2",length(MeanGTR2)),rep("R0=3",length(MeanGTR3)),rep("R0=4",length(MeanGTR4)),rep("R0=5",length(MeanGTR5)) )) 

ggplot(data = dfR0, aes(x=R0, y=MeanRGT, fill=R0))+geom_violin(trim = TRUE)+stat_summary(fun.data = mean_sdl,mult=1,geom = "pointrange",col="black")+scale_fill_brewer(palette="Dark2")+geom_hline(yintercept = 0.5,lty=2,col="black",show.legend = FALSE)+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=11),
  #legend.position = "top",
  legend.title = element_text(size=12, face = "bold"),
  axis.text.x = element_text(size=11),
  axis.title = element_text(size=12),
  axis.text.y = element_text(size=11),
)+ylab(TeX("MeanRGT($\\omega$) - Time Units"))

################################################################
# competition
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R1.5_mu1_IPDistConst_RateConst_.RData")
MeanGTR1.5.c<-comp.overall.eff
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R2_mu1_IPDistConst_RateConst_.RData")
MeanGTR2.c<-comp.overall.eff
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistConst_RateConst_.RData")
MeanGTR3.c<-comp.overall.eff
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R4_mu1_IPDistConst_RateConst_.RData")
MeanGTR4.c<-comp.overall.eff
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R5_mu1_IPDistConst_RateConst_.RData")
MeanGTR5.c<-comp.overall.eff

dfR0.c<-data.frame(MeanRGT.c=c(MeanGTR1.5.c,MeanGTR2.c,MeanGTR3.c,MeanGTR4.c, MeanGTR5.c),R0= c(rep("R0=1.5",length(MeanGTR1.5.c)),rep("R0=2",length(MeanGTR2.c)),rep("R0=3",length(MeanGTR3.c)),rep("R0=4",length(MeanGTR4.c)),rep("R0=5",length(MeanGTR5.c)) )) 
c1<-ggplot(data = dfR0.c, aes(x=R0, y=MeanRGT.c, fill=R0))+geom_violin(trim = TRUE)+stat_summary(fun.data = mean_sdl,mult=1,geom = "pointrange",col="black")+scale_fill_brewer(palette="Greens")+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14),
  axis.title = element_text(size=18),
  axis.text.y = element_text(size=14),
)+ylab(TeX("Competition Intensity ($\\p_{ec}$)"))


# depletion
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R1.5_mu1_IPDistConst_RateConst_.RData")
MeanGTR1.5.d<-abs(depletion.intensity)
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R2_mu1_IPDistConst_RateConst_.RData")
MeanGTR2.d<-abs(depletion.intensity)
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistConst_RateConst_.RData")
MeanGTR3.d<-abs(depletion.intensity)
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R4_mu1_IPDistConst_RateConst_.RData")
MeanGTR4.d<-abs(depletion.intensity)
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R5_mu1_IPDistConst_RateConst_.RData")
MeanGTR5.d<-abs(depletion.intensity)

dfR0.d<-data.frame(MeanRGT.d=c(MeanGTR1.5.d,MeanGTR2.d,MeanGTR3.d,MeanGTR4.d, MeanGTR5.d),R0= c(rep("R0=1.5",length(MeanGTR1.5.d)),rep("R0=2",length(MeanGTR2.d)),rep("R0=3",length(MeanGTR3.d)),rep("R0=4",length(MeanGTR4.d)),rep("R0=5",length(MeanGTR5.d)) )) 
c2<-ggplot(data = dfR0.d, aes(x=R0, y=MeanRGT.d, fill=R0))+geom_violin(trim = TRUE)+stat_summary(fun.data = mean_sdl,mult=1,geom = "pointrange",col="black")+scale_fill_brewer(palette="Blues")+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14),
  axis.title = element_text(size=18),
  axis.text.y = element_text(size=14),
)+ylab(TeX("Depletion Intesity ($\\varphi$)"))


ggarrange(c1,c2, labels = c("A","B"), ncol = 2, nrow = 1) #saved: 12/5

######################################################################
#### IP
######################################################################
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistConst_RateConst_.RData")
MeanGTR1.5<-MeanGT
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistGamma(2,0.5)_RateConst_.RData")
MeanGTR2<-MeanGT
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistExp_RateConst_.RData")
MeanGTR3<-MeanGT
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistGamma(0.5,2)_RateConst_.RData")
MeanGTR4<-MeanGT

dfR0<-data.frame(MeanRGT=c(MeanGTR1.5,MeanGTR2,MeanGTR3,MeanGTR4),IPD= c(rep("Const(1)",length(MeanGTR1.5)),rep("Gamma(2,0.5)",length(MeanGTR2)),rep("Exp(1)",length(MeanGTR3)),rep("Gamma(0.5,2)",length(MeanGTR4)))) 
dfR0$IPD<-factor(dfR0$IPD, levels=c("Const(1)","Gamma(2,0.5)","Exp(1)","Gamma(0.5,2)"))


ggplot(data = dfR0, aes(x=IPD, y=MeanRGT, fill=IPD))+geom_violin(trim = TRUE)+stat_summary(fun.data = mean_sdl,mult=1,geom = "pointrange",col="black")+scale_fill_brewer(palette="Set3")+xlab("Infectious Period Distribution (IPD) - Variance")+scale_x_discrete(labels=c("0","0.5","1","2")) +geom_hline(yintercept = 0.5,lty=2,color="#8DD3C7",lwd=1.1)+geom_hline(yintercept = 0.75,lty=2,color="#FFFFB3",lwd=1.1)+geom_hline(yintercept = 1,lty=2,color="#BEBADA",lwd=1.1)+geom_hline(yintercept = 1.5,lty=2,color="#FB8072",lwd=1.1)+ylim(c(0,1.5))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=11),
  #legend.position = "top",
  legend.title = element_text(size=12, face = "bold"),
  axis.text.x = element_text(size=11),
  axis.title = element_text(size=12),
  axis.text.y = element_text(size=11),
)+ylab(TeX("MeanRGT ($\\omega$) - Time Units"))


################################################################
# competition
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistConst_RateConst_.RData")
MeanGTR1.c<-comp.overall.eff
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistGamma(2,0.5)_RateConst_.RData")
MeanGTR2.c<-comp.overall.eff
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistExp_RateConst_.RData")
MeanGTR3.c<-comp.overall.eff
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistGamma(0.5,2)_RateConst_.RData")
MeanGTR4.c<-comp.overall.eff

dfR0.c<-data.frame(MeanRGT=c(MeanGTR1.c,MeanGTR2.c,MeanGTR3.c,MeanGTR4.c),IPD= c(rep("Const(1)",length(MeanGTR1.c)),rep("Gamma(2,0.5)",length(MeanGTR2.c)),rep("Exp(1)",length(MeanGTR3.c)),rep("Gamma(0.5,2)",length(MeanGTR4.c)))) 
dfR0.c$IPD<-factor(dfR0.c$IPD, levels=c("Const(1)","Gamma(2,0.5)","Exp(1)","Gamma(0.5,2)"))
c1<-ggplot(data = dfR0.c, aes(x=IPD, y=MeanRGT, fill=IPD))+geom_violin(trim = FALSE)+stat_summary(fun.data = mean_sdl,mult=1,geom = "pointrange",col="black")+scale_fill_brewer(palette="Greens")+xlab("Infectious Period Distribution (IPD) - Variance")+scale_x_discrete(labels=c("0","0.5","1","2"))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14),
  axis.title = element_text(size=18),
  axis.text.y = element_text(size=14),
)+ylab(TeX("Competition Intesity ($\\p_{ec}$)"))

plot(c1)

# depletion
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistConst_RateConst_.RData")
MeanGTR1.5.d<-abs(depletion.intensity)
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistExp_RateConst_.RData")
MeanGTR2.d<-abs(depletion.intensity)
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistExp_RateConst_.RData")
MeanGTR3.d<-abs(depletion.intensity)
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistGamma(0.5,2)_RateConst_.RData")
MeanGTR4.d<-abs(depletion.intensity)

dfR0.d<-data.frame(MeanRGT.d=c(MeanGTR1.5.d,MeanGTR2.d,MeanGTR3.d,MeanGTR4.d),IPD= c(rep("Const(1)",length(MeanGTR1.5.d)),rep("Gamma(2,0.5)",length(MeanGTR2.d)),rep("Exp(1)",length(MeanGTR3.d)),rep("Gamma(0.5,2)",length(MeanGTR4.d)))) 
dfR0.d$IPD<-factor(dfR0.d$IPD, levels=c("Const(1)","Gamma(2,0.5)","Exp(1)","Gamma(0.5,2)"))
c2<-ggplot(data = dfR0.d, aes(x=IPD, y=MeanRGT.d, fill=IPD))+geom_violin(trim = FALSE)+stat_summary(fun.data = mean_sdl,mult=1,geom = "pointrange",col="black")+scale_fill_brewer(palette="Blues")+xlab("Infectious Period Distribution (IPD) - Variance")+scale_x_discrete(labels=c("0","0.5","1","2"))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14),
  axis.title = element_text(size=18),
  axis.text.y = element_text(size=14),
)+ylab(TeX("Depletion Intesity ($\\varphi$)"))
plot(c2)

ggarrange(c1,c2, labels = c("A","B"), ncol = 2, nrow = 1) #saved: 12/5

######################################################################
#### Population Size
######################################################################
MeanGTR1.5<-MeanGT
MeanGTR2<-MeanGT
MeanGTR3<-MeanGT
MeanGTR4<-MeanGT

dfR0<-data.frame(MeanRGT=c(MeanGTR1.5,MeanGTR2,MeanGTR3,MeanGTR4),N= c(rep("4",length(MeanGTR1.5)),rep("20",length(MeanGTR2)),rep("100",length(MeanGTR3)),rep("1000",length(MeanGTR4)))) 
dfR0$N<-factor(dfR0$N, levels=c("4","20","100","1000"))
ggplot(data = dfR0, aes(x=N, y=MeanRGT, fill=N))+geom_violin(trim = TRUE)+ theme_minimal()+stat_summary(fun.data = mean_sdl,mult=1,geom = "pointrange",col="red")+scale_fill_brewer(palette="Set2")+xlab("Population size (N)")+geom_hline(yintercept = 0.5, lwd=0.9,lty=2)+ylab(TeX("MeanRGT ($\\omega$) - Time Unit"))

ggplot(data = dfR0, aes(x=N, y=MeanRGT, fill=N))+geom_violin(trim = FALSE)+stat_summary(fun.data = mean_sdl,mult=1,geom = "pointrange",col="red")+scale_fill_brewer(palette="Set2")+xlab("Population size (N)")+scale_x_discrete(labels=c("4","20","100","1000"))+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=11),
  #legend.position = "top",
  legend.title = element_text(size=12, face = "bold"),
  axis.text.x = element_text(size=11),
  axis.title = element_text(size=12),
  axis.text.y = element_text(size=11),
)+ylab(TeX("MeanRGT ($\\omega$) - Time Units"))



################################################################
# competition
MeanGTR1.5.c<-comp.overall.eff
MeanGTR2.c<-comp.overall.eff
MeanGTR3.c<-comp.overall.eff
MeanGTR4.c<-comp.overall.eff

dfR0.c<-data.frame(MeanRGT.c=c(MeanGTR1.5.c,MeanGTR2.c,MeanGTR3.c,MeanGTR4.c),N= c(rep("4",length(MeanGTR1.5.c)),rep("20",length(MeanGTR2.c)),rep("100",length(MeanGTR3.c)),rep("1000",length(MeanGTR4.c)))) 
dfR0.c$N<-factor(dfR0.c$N, levels=c("4","20","100","1000"))
c1<-ggplot(data = dfR0.c, aes(x=N, y=MeanRGT.c, fill=N))+geom_violin(trim = TRUE)+ theme_minimal()+stat_summary(fun.data = mean_sdl,mult=1,geom = "pointrange",col="red")+scale_fill_brewer(palette="Greens")+ylab("Proportion Effective Competition Events")+xlab("Population size (N)")
ggplot(data = dfR0.c, aes(x=N, y=MeanRGT.c, fill=N))+geom_violin(trim = TRUE)+ theme_minimal()+stat_summary(fun.data = mean_sdl,mult=1,geom = "pointrange",col="black")+scale_fill_brewer(palette="Greens")+ylab("Proportion Effective Competition Events")+xlab("Population size (N)")
plot(c1)


# depletion
MeanGTR1.5.d<-abs(depletion.intensity)
MeanGTR2.d<-abs(depletion.intensity)
MeanGTR3.d<-abs(depletion.intensity)
MeanGTR4.d<-abs(depletion.intensity)

dfR0.d<-data.frame(MeanRGT.d=c(MeanGTR1.5.d,MeanGTR2.d,MeanGTR3.d,MeanGTR4.d),N= c(rep("4",length(MeanGTR1.5.d)),rep("20",length(MeanGTR2.d)),rep("100",length(MeanGTR3.d)),rep("1000",length(MeanGTR4.d)))) 
dfR0.d$N<-factor(dfR0$N, levels=c("4","20","100","1000"))
c2<-ggplot(data = dfR0.d, aes(x=N, y=MeanRGT.d, fill=N))+geom_violin(trim = TRUE)+ theme_minimal()+stat_summary(fun.data = mean_sdl,mult=1,geom = "pointrange",col="red")+scale_fill_brewer(palette="Blues")+ylab("Depletion Intensity")+xlab("Population size (N)")
ggplot(data = dfR0.d, aes(x=N, y=MeanRGT.d, fill=N))+geom_violin(trim = TRUE)+ theme_minimal()+scale_fill_brewer(palette="Blues")+ylab("Depletion Intensity")+xlab("Population size (N)")
ggarrange(c1,c2, labels = c("A","B"), ncol = 2, nrow = 1) #saved: 12/5


######################################################################
#### Rates
######################################################################
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/Rates/CompDepl_n1000_R3_mu1_IPDistConst_RateInf_.RData")
MeanGTR1.5<-MeanGT
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/Rates/CompDepl_n1000_R3_mu1_IPDistConst_RateBC_.RData")
MeanGTR2<-MeanGT
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistConst_RateConst_.RData")
MeanGTR3<-MeanGT

dfR0<-data.frame(MeanRGT=c(MeanGTR1.5,MeanGTR2,MeanGTR3),Rate= c(rep("Peaked",length(MeanGTR1.5)),rep("Decreasing",length(MeanGTR2)),rep("Constant",length(MeanGTR3)))) 
dfR0$Rate<-factor(dfR0$Rate, levels=c("Peaked","Decreasing","Constant"))
ggplot(data = dfR0, aes(x=Rate, y=MeanRGT, fill=Rate))+geom_violin(trim = FALSE)+stat_summary(fun.data = mean_sdl,mult=1,geom = "pointrange",col="black")+scale_fill_brewer(palette="Set1")+xlab("Effective Contact Profile")+geom_hline(yintercept = 0.5, lty=2, lwd=0.9, color="#4DAF4A")+geom_hline(yintercept = 0.21, lty=2, lwd=0.9, color="#377EB8")+geom_hline(yintercept = 0.502, lty=2, lwd=0.9, color="#E41A1C")+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=12),
  #legend.position = "top",
  legend.title = element_text(size=12, face = "bold"),
  axis.text.x = element_text(size=11),
  axis.title = element_text(size=12),
  axis.text.y = element_text(size=11),
)+ylab(TeX("MeanRGT ($\\omega$) - Time Units"))






################################################################
# competition
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/Rates/CompDepl_n1000_R3_mu1_IPDistConst_RateInf_.RData")
MeanGTR1.5.c<-comp.overall.eff
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/Rates/CompDepl_n1000_R3_mu1_IPDistConst_RateBC_.RData")
MeanGTR2.c<-comp.overall.eff
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistConst_RateConst_.RData")
MeanGTR3.c<-comp.overall.eff

dfR0.c<-data.frame(MeanRGT.c=c(MeanGTR1.5.c,MeanGTR2.c,MeanGTR3.c),Rate= c(rep("Peaked",length(MeanGTR1.5.c)),rep("Decreasing",length(MeanGTR2.c)),rep("Constant",length(MeanGTR3.c)))) 
dfR0.c$Rate<-factor(dfR0.c$Rate, levels=c("Peaked","Decreasing","Constant"))
c1<-ggplot(data = dfR0.c, aes(x=Rate, y=MeanRGT.c, fill=Rate))+geom_violin(trim = TRUE)+stat_summary(fun.data = mean_sdl,mult=1,geom = "pointrange",col="black")+scale_fill_brewer(palette="Greens")+ylab("Proportion Effective Competition Events")+xlab("Effective Contact Profile")+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14),
  axis.title = element_text(size=18),
  axis.text.y = element_text(size=14),
)+ylab(TeX("Competition Intesity ($\\p_{ec}$)"))
plot(c1)


# depletion
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/Rates/CompDepl_n1000_R3_mu1_IPDistConst_RateInf_.RData")
MeanGTR1.5.d<-abs(depletion.intensity)
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/Rates/CompDepl_n1000_R3_mu1_IPDistConst_RateBC_.RData")
MeanGTR2.d<-abs(depletion.intensity)
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistConst_RateConst_.RData")
MeanGTR3.d<-abs(depletion.intensity)

dfR0.d<-data.frame(MeanRGT.d=c(MeanGTR1.5.d,MeanGTR2.d,MeanGTR3.d),Rate= c(rep("Peaked",length(MeanGTR1.5.d)),rep("Decreasing",length(MeanGTR2.d)),rep("Constant",length(MeanGTR3.d)))) 
dfR0.d$Rate<-factor(dfR0.d$Rate, levels=c("Peaked","Decreasing","Constant"))
c2<-ggplot(data = dfR0.d, aes(x=Rate, y=MeanRGT.d, fill=Rate))+geom_violin(trim = TRUE)+ stat_summary(fun.data = mean_sdl,mult=1,geom = "pointrange",col="black")+scale_fill_brewer(palette="Blues")+ylab("Depletion Intensity")+xlab("Effective Contact Profile")+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14),
  axis.title = element_text(size=18),
  axis.text.y = element_text(size=14),
)+ylab(TeX("Depletion Intesity ($\\varphi$)"))
plot(c2)
ggarrange(c1,c2, labels = c("A","B"), ncol = 2, nrow = 1) #saved: 12/5



#####################################################################
# Competition
#####################################################################

MeanGTR1.5<-MeanGT[which(!is.na(MeanGT))]
MeanGTR2<-MeanGT[which(!is.na(MeanGT))]
MeanGTR3<-MeanGT[which(!is.na(MeanGT))]

dfR0<-data.frame(MeanRGT=c(MeanGTR1.5,MeanGTR2,MeanGTR3),Setting= c(rep("A",length(MeanGTR1.5)),rep("B",length(MeanGTR2)),rep("C",length(MeanGTR3)))) 
dfR0$Setting<-factor(dfR0$Setting, levels=c("A","B","C"))
ggplot(data = dfR0, aes(x=Setting, y=MeanRGT, fill=Setting))+geom_violin(trim = TRUE)+ theme_minimal()+stat_summary(fun.data = mean_sdl,mult=1,geom = "pointrange",col="black")+scale_fill_brewer(palette="Set2")+xlab("Index Cases")+scale_x_discrete(labels=c("0","0.5","1","2"))



################################################################3
# SUPPLEMENTARY MATERIAL
################################################################3
# 
# FWD AND BKWD GEN TIME
#
n<-1000 #pop size
mu<-1 #length of the infectious period. It is assumed to be constant
#lambda<-1 #daily contact rate
R0<-4 # reproduction number
#seed<-425513
source("SimFunctionReGTs.R")

#CONTACT-RATE
#run the simulations
epi<-list()
set.seed(62134)
nSim<-10
for (i in 1:nSim){
  epi[[i]]<-sim.comp.depl.hom.k(mu=mu,lambda = R0, n=n)
  print(i)
}

#simulations that are not extinct 
not.extinct<-0
extinct<-max(1,floor(0.1*(n)))
for (i in 1:nSim){
  ifelse(length(epi[[i]]$generation$CurrentTime)>extinct,not.extinct[i]<-i,not.extinct[i]<-NA)
}

################################################3##############3
#select one simulation
w<-7

epi.sim<-epi[[w]]
#Competition events and effects
# in a data frame is reported for every generations affected by competition whether is effective or not, the mean among the contacts proposed by the potential infector, the deviation from this mean value and the number of infectors
comp.eff<-data.frame("Time.of.generation"=0,"is.effective"=0,"mean.proposed.contact"=0,"potential.gain"=0, "Potential.infectors"=0)
time.of.gen<-0

for (i in 1:n){
  time.of.gen[i]<-epi.sim$competition[[i]][1]
}


current.gen<-0
is.eff<-0
infs<-0
m<-0
expgt<-0.5
pg<-0
infectee.generation<-epi.sim$time.events[which(epi.sim$time.events[,2]==1),]
comp.eff[1,]<-c(0,0,0,0,0)
for (i in 2:length(infectee.generation[,1])){
  current.gen<-infectee.generation[i,3]
  ifelse(epi.sim$competition[[current.gen]][2]==min(epi.sim$competition[[current.gen]][2:length(epi.sim$competition[[current.gen]])]) & length(epi.sim$competition[[current.gen]])>2,is.eff<-1,is.eff<-0)
  m<-mean(epi.sim$competition[[current.gen]][2:length(epi.sim$competition[[current.gen]])])
  #pg<-(m-epi.sim$competition[[current.gen]][2])/m
  
  pg<-(expgt-epi.sim$competition[[current.gen]][2])/expgt
  ifelse(length(epi.sim$competition[[current.gen]])>2, infs<-(length(epi.sim$competition[[current.gen]])-1),infs<-NA)
  comp.eff[i,]<-c(epi.sim$competition[[current.gen]][1],is.eff,m,pg,infs)
}
generation.aff.competition<-length(which(comp.eff$is.effective==1))/length(comp.eff$is.effective)
gene.comp.eff<-which(comp.eff$is.effective==1)
total.mean.gain<-mean(comp.eff$potential.gain[gene.comp.eff])
mean.infectors<-mean(comp.eff$Potential.infectors, na.rm = T)

#How many competition events in a small time interval
tEnd<-epi.sim$time.events[length(epi.sim$time.events[,1])]
t<-seq(0,tEnd,0.05)
NCompEvnts<-NULL
tc<-comp.eff[which(comp.eff$is.effective==1),]
for (i in t){
  f<-intersect(which(tc$Time.of.generation < i+0.1), which(tc$Time.of.generation >i))
  NCompEvnts<-c(NCompEvnts, length(f))
}

#Depletion intensity
library("drc")
max.depletion<-NULL
epi.curve.1<-0
suscept<-0
#compute the number of susceptibles over time
time.events<-epi.sim$time.events
for (s in 1:length(time.events[,1]) ){
  epi.curve.1[s]<-length(which(time.events[1:s,2]==1))-length(which(time.events[1:s,2]==0.1))
  suscept[s]<-n-epi.curve.1[s]-length(which(time.events[1:s,2]==0.1))
}
suscept<-suscept/n # I just look at the proportion, i.e. the probability of encountering a susceptible given a contact
mL <- drm(suscept~time.events[,1], fct = L.5(), type = "continuous") #Fit a logistic curve
b<-as.numeric(mL$coefficients[1])
c<-as.numeric(mL$coefficients[2])
d<-as.numeric(mL$coefficients[3])
e<-as.numeric(mL$coefficients[4])
f<-as.numeric(mL$coefficients[5])

b.vec<-rep(b,length(time.events[,1]))
c.vec<-rep(c,length(time.events[,1]))
d.vec<-rep(d,length(time.events[,1]))
e.vec<-rep(e,length(time.events[,1]))
f.vec<-rep(f,length(time.events[,1]))
unit<-rep(1,length(time.events[,1]))

susc.der<--(d.vec-c.vec)*(unit+exp(b.vec*(time.events[,1]-e.vec)))^(-f.vec-unit)*b.vec*exp(b.vec*(time.events[,1]-e.vec)) #derivative

susc.der1<-function(t){
  return(-(d-c)*(1+exp(b*(t-e)))^(-f-1)*b*exp(b*(t-e)))
}

max.depletion<-c(time.events[which(susc.der==min(susc.der)),1], min(susc.der))

#epidemic curve (proportion)
int<-seq(0,tEnd,0.05)
time.events<-epi.sim$time.events
infection<-time.events[which(time.events[,2]==1.1),]
epi.curve.1<-0
suscept<-0
new.inf<-NULL
for (i in 1:length(time.events[,1]) ){
  epi.curve.1[i]<-length(which(time.events[1:i,2]==1))-length(which(time.events[1:i,2]==0.1))
  suscept[i]<-n-epi.curve.1[i]-length(which(time.events[1:i,2]==0.1))
}
for (i in int){
  te<-intersect(which(time.events[,1]<(i+0.03)), which(time.events[,1]>i))
  new.inf<-c(new.inf,length(which(time.events[te,2]==1)))
}

individual.gen<-matrix(NA,nrow = length(unique(epi.sim$generation[,2])),ncol = 2) # first colum infectious time and second column mean among the "succesfull" effective contacts
for (i in 1:length(unique(epi.sim$generation[,2]))){
  infector<-unique(epi.sim$generation[,2])[i]
  individual.gen[i,1]<-epi.sim$status.matrix[infector,2]
  individual.gen[i,2]<-mean(epi.sim$generation$EffectiveContact[which(epi.sim$generation$Infector==infector)])
}
temp<-individual.gen #reorder the data
temp1<-sort(temp[,1])
for (i in 1:length(individual.gen[,1])) {
  individual.gen[i,1]<-temp1[i]
  individual.gen[i,2]<-temp[which(temp1[i]==temp[,1]),2]
}


##########################################################
##########################################################
# Fig S1 - Forward Gen time
fwd.gen<-fwd.gen.int(epi.sim$status.matrix,epi.sim$time.events)
fwd.gen<-fwd.gen[!is.na(fwd.gen[,2]),]
spline.fwd<-loess.smooth(fwd.gen[,1],fwd.gen[,2])
EffKFWD<-loess.smooth(individual.gen[,1], individual.gen[,2])
k11<-mu/(1+R0)


par(mar=c(3,5,2,3))

plot(time.events[,1], suscept/n, ylim = c(0,1),xlim = c(0,tEnd), type = "l", ylab = "Time Units ", xlab = "Days", lwd=2)
abline(v=c(max.depletion[1]+k11,max.depletion[1]+2*k11,max.depletion[1]+3*k11,max.depletion[1]+4*k11), lty=2, col="grey")
abline(v=max.depletion[1], col="purple", lty=2, lwd=2)
par(new=T)
plot(EffKFWD$x,k11*EffKFWD$y, col="blue", xlim = c(0,tEnd), ylim = c(0,1), type = "l", xlab = " ", xaxt="n", ylab = " ", yaxt="n", lwd=2)

plot(EffKFWD$x,k11*EffKFWD$y, col="blue", xlim = c(0,3.5), ylim = c(0.1,0.8), type = "l", xlab = "Time Units", ylab = "Time Units", lwd=2, cex.axis=1.2, cex.lab=1.2)

par(new=T)
par(new=T)
plot(spline.fwd$x,spline.fwd$y, ylim = c(0,1), xlim=c(0,tEnd),col="darkgreen",type = "l", xlab = " ", xaxt= "n", ylab = "", yaxt="n", lty=1, lwd=2)
axis(4)
mtext('Proportion Susceptibles',4,line=2)
#par(new=T)
#plot(x=comp.eff$Time.of.generation[which(comp.eff$is.effective==1)], y=rep(0.1,length(comp.eff$Time.of.generation[which(comp.eff$is.effective==1)]) ), xlim = c(0, tEnd), ylim = c(0,1))
rbPal <- colorRampPalette(c("yellow","red"))
z <- rbPal(10)[as.numeric(cut(NCompEvnts[which(NCompEvnts>1)],breaks = 10))]
points(t[which(NCompEvnts>1)] , rep(0.1, length(which(NCompEvnts>1))),pch = 20,col = z,lwd=2)
legend("topright", legend = c("GenTime","S","EAG","Comp","MaxDepl"), col=c("darkgreen", "black","blue","red","purple"), lty=c(1,1,1,NA,2),pch = c(NA,NA,NA,20,NA), cex = 0.8)
legend("topright", legend = c("MeanGenTime"), col=c("blue"), lty=c(1), cex = 1.2)
#########################
#Fig S2
bkw.gen<-bkw.gen.int(epi.sim$status.matrix,epi.sim$time.events)[-1,]
smooth.succ.effK<-loess.smooth(epi.sim$generation[-1,1],epi.sim$generation[-1,4], span = 0.75)
smooth.bkw<-loess.smooth(bkw.gen[-1,1],bkw.gen[-1,2], span=0.75) #spline smoother
plot(smooth.succ.effK$x, k11 *smooth.succ.effK$y, main="Backward scheme",col="blue", xlim = c(0,tEnd), type = "l", ylim = c(0,1), ylab = "Days", xlab = "Days",lwd=2)

plot(int, new.inf, ylim = c(0,50),xlim = c(0,tEnd), type = "l", ylab = " ", xlab = "Days", lwd=1, yaxt="n")
#plot(time.events[,1], epi.curve.1/n, ylim = c(0,1),xlim = c(0,tEnd), type = "l", ylab = " ", xlab = "Days", lwd=2)
abline(v=max.depletion[1], col="purple", lty=2, lwd=2)
par(new=T)
plot(smooth.succ.effK$x, k11 *smooth.succ.effK$y, main=" ",col="blue", xlim = c(0,tEnd), type = "l", ylim = c(0,1), ylab = " Time Units", xlab = "Days",lwd=2)
par(new=T)
plot(smooth.bkw$x,smooth.bkw$y, ylim = c(0,1), xlim=c(0,tEnd),col="darkgreen",type = "l", xlab = " ", xaxt="n", ylab = " ", yaxt="n",lty=1,lwd=2)
#par(new=T)
axis(4)
mtext('Proportion Infectives',4,line=2)

#plot(x=comp.eff$Time.of.generation[which(comp.eff$is.effective==1)], y=rep(0.1,length(comp.eff$Time.of.generation[which(comp.eff$is.effective==1)]) ), xlim = c(0, tEnd), ylim = c(0,1))
rbPal <- colorRampPalette(c("yellow","red"))
z <- rbPal(10)[as.numeric(cut(NCompEvnts[which(NCompEvnts>1)],breaks = 10))]
points(t[which(NCompEvnts>1)] , rep(0.1, length(which(NCompEvnts>1))),pch = 20,col = z,lwd=2, xlim=c(0,tEnd))
legend("topright", legend = c("GenTime","NewCases","EAG","Comp","MaxDepl"), col=c("darkgreen", "black","blue","red","purple"), lty=c(1,1,1,NA,2),pch = c(NA,NA,NA,20,NA), cex = 0.8)


#################################################3
#Correlation plots
#################################################3


data.cor<-data.frame("Competition"=comp.overall.eff,"Depletion"=abs(depletion.intensity))

cR15<-comp.overall.eff
cR2<-comp.overall.eff
cR3<-comp.overall.eff
cR4<-comp.overall.eff
cR5<-comp.overall.eff

dR15<-abs(depletion.intensity)
dR2<-abs(depletion.intensity)
dR3<-abs(depletion.intensity)
dR4<-abs(depletion.intensity)
dR5<-abs(depletion.intensity)


kol<-brewer.pal(n = 5, name = "Dark2")

library(latex2exp)


#R15
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R1.5_mu1_IPDistConst_RateConst_.RData")
data.cor<-data.frame("Competition"=comp.overall.eff,"Depletion"=abs(depletion.intensity))
a1<-ggplot(data.cor) + aes(x = Competition, y = Depletion) + geom_point(colour = kol[1]) +xlab(TeX("$p_{ec}$"))+ylab(TeX("$\\varphi$"))+theme(
  panel.background = element_rect(fill = "white", colour = "white",size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.05, linetype = 'solid',colour = "gray"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "white"),
  axis.text.x = element_text(size=16),
  axis.title = element_text(size=20),
  axis.text.y = element_text(size=16),
)+scale_y_continuous(breaks = c(0,0.25,0.5,0.75))
val1<- cor(data.cor$Competition,data.cor$Depletion, method = "spearman")

#R2
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R2_mu1_IPDistConst_RateConst_.RData")
data.cor<-data.frame("Competition"=comp.overall.eff,"Depletion"=abs(depletion.intensity))
a2<-ggplot(data.cor) + aes(x = Competition, y = Depletion) + geom_point(colour = kol[2]) +xlab(TeX("$p_{ec}$"))+ylab(TeX("$\\varphi$"))+theme(
  panel.background = element_rect(fill = "white", colour = "white",size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.05, linetype = 'solid',colour = "gray"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "white"),
  axis.text.x = element_text(size=16),
  axis.title = element_text(size=20),
  axis.text.y = element_text(size=16),
)+scale_y_continuous(breaks = c(0.2,0.4,0.6))+scale_x_continuous(breaks = c(0.075,0.125,0.1,0.15,0.175))
val2<- cor(data.cor$Competition,data.cor$Depletion, method = "spearman")

#R3
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R3_mu1_IPDistConst_RateConst_.RData")
data.cor<-data.frame("Competition"=comp.overall.eff,"Depletion"=abs(depletion.intensity))
a3<-ggplot(data.cor) + aes(x = Competition, y = Depletion) + geom_point(colour = kol[3]) +xlab(TeX("$p_{ec}$"))+ylab(TeX("$\\varphi$"))+theme(
  panel.background = element_rect(fill = "white", colour = "white",size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.05, linetype = 'solid',colour = "gray"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "white"),
  axis.text.x = element_text(size=16),
  axis.title = element_text(size=20),
  axis.text.y = element_text(size=16),
)
val3<- cor(data.cor$Competition,data.cor$Depletion, method = "spearman")

#R4
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R4_mu1_IPDistConst_RateConst_.RData")
data.cor<-data.frame("Competition"=comp.overall.eff,"Depletion"=abs(depletion.intensity))
a4<-ggplot(data.cor) + aes(x = Competition, y = Depletion) + geom_point(colour = kol[4]) +xlab(TeX("$p_{ec}$"))+ylab(TeX("$\\varphi$"))+theme(
  panel.background = element_rect(fill = "white", colour = "white",size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.05, linetype = 'solid',colour = "gray"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "white"),
  axis.text.x = element_text(size=16),
  axis.title = element_text(size=20),
  axis.text.y = element_text(size=16),
)+scale_y_continuous(breaks = c(0.5,1,1.5))
val4<- cor(data.cor$Competition,data.cor$Depletion, method = "spearman")

#R5
load("~/Documents/Work/PhD/RealizeGenTime/NewReGts/SimSummary/R0/CompDepl_n1000_R5_mu1_IPDistConst_RateConst_.RData")
data.cor<-data.frame("Competition"=comp.overall.eff,"Depletion"=abs(depletion.intensity))
a5<-ggplot(data.cor) + aes(x = Competition, y = Depletion) + geom_point(colour = kol[5]) +xlab(TeX("$p_{ec}$"))+ylab(TeX("$\\varphi$"))+theme(
  panel.background = element_rect(fill = "white", colour = "white",size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.05, linetype = 'solid',colour = "gray"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "white"),
  axis.text.x = element_text(size=16),
  axis.title = element_text(size=20),
  axis.text.y = element_text(size=16),
)+scale_y_continuous(breaks = c(0.5,1,1.5))
val5<- cor(data.cor$Competition,data.cor$Depletion, method = "spearman")

ggarrange(a1,a2,a3,a4,a5, labels = c("R=1.5","R=2","R=3","R=4","R=5"),  font.label = list(size = 20, color = "black", face = "bold", family = NULL),ncol = 2, nrow = 3) #saved: 12/5









