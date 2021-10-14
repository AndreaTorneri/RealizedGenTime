args <- commandArgs(trailingOnly = TRUE)

out = args[1] #working directory
cat(",out=",out)
cores = as.numeric(args[2]) #number of cores to run in parallel
cat(",cores=",cores)
R = as.numeric(args[3]) #basic Reproduction number
cat(",R=",R)
n = as.numeric(args[4]) # population size
cat(",n=",n)
mu = as.numeric(args[5]) # infectious period length(average)
cat(",mu=",mu)
n.InInf = as.numeric(args[6]) # number of initial infectious individuals
cat(",n.InInf=",n.InInf)
n.Rec = as.numeric(args[7]) # number of initial recovered
cat(",n.Rec=",n.Rec)
tv = as.numeric(args[8]) # time of vaccination
cat(",tv=",tv)
p = as.numeric(args[9]) # vaccination coverage
cat(",p=",p)


print(paste("R",R,"_n",n,"_mu",mu,"_nInf",n.InInf,"_nR",n.Rec,"_tv",tv,"_p",p, sep = ""))


IPDIST<-"Const"
RATE<-"Const"

source("SimFunctionReGTs.R")



#lambda.g<-8.29
#SettingInfectivityParameters
# direct computation
#source("R_comp_netw.R")
#ratio_hhgl<-1.31/lambda.g*lambda.h
#R.rif<-R.1
#nSim<-100
#tol<-0.05
#trs.prms<-R0.comp(ratio_hhgl=ratio_hhgl, HH.network = HH.network, nSim = nSim, tol=tol,R.rif = R.rif )






# Parameters
#library("drc")


#For the competition scenario
nSim<-3000
epi.outbreak<-list()
nSeed<-14022020
set.seed(nSeed)

for (i in 1:nSim){
   print(i)
  #epi.outbreak[[i]]<-sim.comp.depl.hom.k(mu=mu,lambda = R, n=n)
  #epi.outbreak[[i]]<-sim.comp.depl.hom.inf(mu=mu,lambda = R, n=n)
  #epi.outbreak[[i]]<-sim.comp.depl.hom.bc(mu=mu,lambda = R, n=n)
  #epi.outbreak[[i]]<-generation.interval.parallel.fixedpropsuscept.hhscen(R=R,n=n,mu=mu,n.InInf = n.InInf,n.Rec = n.Rec)
  epi.outbreak[[i]]<-sim.comp.depl.vacc(mu=mu,lambda = R,n=n,tv=tv,p=p)
  }


#finalSize<-NULL
#not.extinct<-NULL
#for (i in 1:nSim){
#  finalSize[i]<-length(which(epi.outbreak[[i]]$time.events[,2]==1.0))
#     if (finalSize[i]>(round(n*0.1))){not.extinct<-c(not.extinct,i)}
  #if (finalSize[i]>1){not.extinct<-c(not.extinct,i)}
#}
#FinSize<-finalSize[not.extinct]

#MeanGT<-NULL
#MeanT<-NULL
#for (i in not.extinct){
#  MeanGT<-c(MeanGT, mean(epi.outbreak[[i]]$GT))
#  MeanT<-c(MeanT, max(epi.outbreak[[i]]$time.events[,1]))
#}

# max.depletion<-matrix(NA, nrow = 1, ncol = 2)
# for (j in not.extinct){
#   epi.curve.1<-0
#   suscept<-0
#   #compute the number of susceptibles over time
#   time.events<-epi.outbreak[[j]]$time.events
#   for (s in 1:length(time.events[,1]) ){
#     epi.curve.1[s]<-length(which(time.events[1:s,2]==1))-length(which(time.events[1:s,2]==0.1))
#     suscept[s]<-n-epi.curve.1[s]-length(which(time.events[1:s,2]==0.1))
#   }
#   suscept<-suscept/n
#   mL <- drm(suscept~time.events[,1], fct = L.5(), type = "continuous") #Fit a logistic curve
#   b<-as.numeric(mL$coefficients[1])
#   c<-as.numeric(mL$coefficients[2])
#   d<-as.numeric(mL$coefficients[3])
#   e<-as.numeric(mL$coefficients[4])
#   f<-as.numeric(mL$coefficients[5])
#   b.vec<-rep(b,length(time.events[,1]))
#   c.vec<-rep(c,length(time.events[,1]))
#   d.vec<-rep(d,length(time.events[,1]))
#   e.vec<-rep(e,length(time.events[,1]))
#   f.vec<-rep(f,length(time.events[,1]))
#   unit<-rep(1,length(time.events[,1]))
##   
#   susc.der<--(d.vec-c.vec)*(unit+exp(b.vec*(time.events[,1]-e.vec)))^(-f.vec-unit)*b.vec*exp(b.vec*(time.events[,1]-e.vec))
#   
#   max.depletion<-rbind(max.depletion,c(which(susc.der==min(susc.der[which(is.finite(susc.der))])), min(susc.der[which(is.finite(susc.der))]) ))
#   #o<-o+1
#   #print(j)
# }


#max.depletion<-max.depletion[-1,]


#comp.overall<-NULL
#comp.overall.eff<-NULL

# for (o in not.extinct){
#   comp.eff<-data.frame("Time.of.generation"=0,"is.effective"=0,"mean.proposed.contact"=0,"potential.gain"=0, "Potential.infectors"=0)
#   time.of.gen<-0
#   current.gen<-0
#   is.eff<-0
#   infs<-0
#   m<-0
#   pg<-0
#   infectee.generation<-epi.outbreak[[o]]$time.events[which(epi.outbreak[[o]]$time.events[,2]==1),]
#   comp.eff[1,]<-c(0,0,0,0,0)
#   for (i in 2:length(infectee.generation[,1])){
#     current.gen<-infectee.generation[i,3]
#     ifelse(epi.outbreak[[o]]$competition[[current.gen]][2]==min(epi.outbreak[[o]]$competition[[current.gen]][2:length(epi.outbreak[[o]]$competition[[current.gen]])]) & length(epi.outbreak[[o]]$competition[[current.gen]])>2,is.eff<-1,is.eff<-0)
#     ifelse(length(epi.outbreak[[o]]$competition[[current.gen]])>2, infs<-(length(epi.outbreak[[o]]$competition[[current.gen]])-1),infs<-NA)
#     comp.eff[i,]<-c(epi.outbreak[[o]]$competition[[current.gen]][1],is.eff,m,pg,infs)
#   }
   
#   generation.aff.competition.eff<-length(which(comp.eff$is.effective==1))/finalSize[o]
#   if (generation.aff.competition.eff==0){
#     total.mean.gain<-NA
#     mean.infectors<-NA
#   }
#   comp.overall.eff<-c(comp.overall.eff,generation.aff.competition.eff)
#   }


name<-paste("CompScenario_","n",n,"_R",R,"_mu",mu, "_IPDist",IPDIST, "_Rate",RATE, "_NinF",n.InInf ,"_tv",tv,"_p",p , "_.RData", sep = "")
setwd(out)
save(epi.outbreak, file = name)

