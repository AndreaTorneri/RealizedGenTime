acc.function.inf<-function(t){
  a<-3.797507
  b<-1.611483
  c<-0.6999198
  #int.v<-13.4638
  int.v<-1.15*5
  value<-(a*(t*12)^b*exp(-c*(t*12)))/(int.v)
  return(value)
}

acc.function.inf.2<-function(t){
  #int.v<-1.15*5
  value<-dgamma(t,shape = 3, scale = 0.2)/(5*0.875348)
  return(value)
  
}

acc.function.inf.3<-function(t){
  #int.v<-1.15*5
  value<-dweibull(t,shape = 3.5, scale = 0.56)/(5*0.9995042)
  return(value)
}

acc.function.k<-function(t){
  d<-1
  value<-1/d
  return(value)
}

acc.function.k.rate<-function(t){
  d<-5
  value<-1/d
  return(value)
}

acc.function.bc<-function(t){
  a<-exp(-0.5*(9*t))
  d<-2*(1-exp(-9))
  int.v<-5*0.1098903
  return((a/d)/int.v)
}


sim.comp.depl.hom.k<-function(mu,lambda,n){
  
  status.matrix <- matrix(NA,nrow = n,ncol = 3) 
  generation<-data.frame("CurrentTime"=0, "Infector"=0, "GenerationLength"=0,"EffectiveContact"=0)
  effective.contacts<-list() #keep track of the effective contacts each infectious individual makes
  generation.individual<-list() # keep track of which effective contacts proposed by the specific individual lead to generations
  competition<-list() #keep track of all the effective contact intervals potential infectors propose to a specific susceptible
  for (i in 1:n){
    effective.contacts[[i]]<-0
    generation.individual[[i]]<-NA
    competition[[i]]<-0
  }
  col.name<-c("infected","time.of.infection","infector")
  dimnames(status.matrix)<-list(NULL,col.name)
  status.matrix[,1]<-0 #all the population is initially susceptible
  
  recovery.vector<-rep(NA,n) #vector giving the recovery times
  infectives<-rep(0,n)
  current.time<-0
  rcv<-rep(Inf,n)
  
  # first infected randomly chosen in the population
  first<-sample(1:n, 1)
  generation$Infector[1]<-first
  status.matrix[first,1] <- 1
  status.matrix[first,2] <- 0
  recovery.vector[first]<-current.time+mu
  rcv[first]<-recovery.vector[first]
  infectives[first]<-1
  time.events<-matrix(NA,1,3)
  time.events[1,]<-c(current.time,1,first)
  
  contact.time<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n),"pr.infectee"=rep(NA,n))   #matrix containing the proposed time of the next possible infectious contact (first colum) 
  # with a randomly selected individual (second column)
  proposed.individual<-0
  index.contact<-rep(0,n) # index that selects the individuals for whom a proposed contact is drawn.
  index.contact[first]<-1
  int.time<-0
  temp.contact.time<-0
  T_g<-0
  recovered<-0
  GT<-NULL
  
  #When only the first pathogen is present
  while(sum(infectives, na.rm = TRUE) > 0){ #while there are still infectives
    #Phase 1: individuals that has to, propose a new social contact
    for (i in which(index.contact==1)){ # for all the individuals that has to propose a global contact
      temp.contact.time<-rexp(1,lambda)+current.time# I generate the next interarrival time for individual i
      index.contact[i]<-0
      if (temp.contact.time<recovery.vector[i]){
        contact.time$pr.infectee[i] <-sample(setdiff(1:n,i),1)
        contact.time$pr.ctc[i]<-temp.contact.time # If it is an infectious contact I keep track of the individual and of the time
      }
    }
    
    #Phase 2: identify the next event: possible infection, recovery or the start of the new pathogen infection
    ifelse(length(which(is.na(contact.time$pr.ctc)==FALSE))>0,T_g<-min(contact.time$pr.ctc, na.rm = T),T_g<-Inf)
    R_a<-min(recovery.vector, na.rm = T)
    
    # Phase 2.2 - next event is a proposed infection with pathogen 1
    if (T_g<R_a){
      current.time<-T_g
      infector<-which(contact.time$pr.ctc ==T_g)
      infectee<-sample(setdiff(1:n,infector),1)
      acceptance.rate<-acc.function.k(t=current.time-status.matrix[infector,2])
      is.eff<-0 #we keep track of a contact that is effective
      if (runif(1)<acceptance.rate){ #the contact is effective
        effective.contacts[[infector]]<-c(effective.contacts[[infector]],current.time)
        is.eff<-1
      }
      
      #competition check
      if (is.eff==1 & ((is.na(status.matrix[infectee,1]))| (!(is.na(status.matrix[infectee,1])) & (status.matrix[infectee,1]==1)))){ #if the contact is effective but made with a recovered or a currently infectious individual
        #      if ((status.matrix[infectee,2]-status.matrix[infector,2])<mu & (status.matrix[infector,2]<status.matrix[infectee,2])){ #we check whether the selected infector could be capable of infecting the selected infectee, i.e. if it was infectious when the infectee susceptible
        if ((status.matrix[infectee,2]< rcv[infector]) & (status.matrix[infector,2]<status.matrix[infectee,2])){ #we check whether the selected infector could be capable of infecting the selected infectee, i.e. if it was infectious when the infectee susceptible          
          competition[[infectee]]<-c(competition[[infectee]],(current.time-status.matrix[infector,2]))            
        }
      }
      
      if (!(is.na(status.matrix[infectee,1])) & (status.matrix[infectee,1]==0) & (is.eff==1)){ #if the contact is effective and the contacted individual is susceptible there is transmission
        status.matrix[infectee,1]<-1
        infectives[infectee]<-1
        status.matrix[infectee,2]<-current.time
        status.matrix[infectee,3]<-infector
        recovery.vector[infectee]<-current.time+mu
        rcv[infectee]<-recovery.vector[infectee]
        index.contact[infectee]<-1
        index.contact[infector]<-1
        GT<-c(GT,(current.time-status.matrix[infector,2]))
        effective.contacts[[infectee]][1]<-current.time #the first element of the effective contact list for a specific individual is his/her infection time
        time.events<-rbind(time.events,c(current.time,1,infectee))
        competition[[infectee]]<-c(current.time, (current.time-status.matrix[infector,2]) )# the first element of the competition list is the generation time that infect the individual
        generation.individual[[infector]]<-c(generation.individual[[infector]],(which(effective.contacts[[infector]]==current.time)-1)) #the effective contact that lead to generation (-1 because the first element of EK is the time of infection)
        generation<-rbind(generation,c(current.time,infector,(current.time-status.matrix[infector,2]),(which(effective.contacts[[infector]]==current.time)-1)))
      }else{
        index.contact[infector]<-1
      }
      contact.time$pr.ctc[infector]<-NA
      #Phase 2.3 a recovery occurs  
    }else{
      current.time<-R_a
      recovered<-which(recovery.vector==R_a)
      recovery.vector[recovered]<-NA
      status.matrix[recovered,1]<-NA
      time.events<-rbind(time.events,c(current.time,0.1,recovered))
      contact.time[recovered,2:5]<-rep(NA,4)
      infectives[recovered]<-NA
      effective.contacts[[recovered]]<-c(effective.contacts[[recovered]],current.time) # the last value of effective contact is the day of recover
    }
  }
  timev.name<-c("time","event","who")
  dimnames(time.events)<-list(NULL,timev.name)
  
  return(list(time.events=time.events, status.matrix=status.matrix, generation=generation, effective.contacts=effective.contacts, generation.individual=generation.individual, competition=competition, GT=GT))
}

sim.comp.depl.hom.inf<-function(mu,lambda,n){
  lambda<-lambda*5
  status.matrix <- matrix(NA,nrow = n,ncol = 3) 
  generation<-data.frame("CurrentTime"=0, "Infector"=0, "GenerationLength"=0,"EffectiveContact"=0)
  effective.contacts<-list() #keep track of the effective contacts each infectious individual makes
  generation.individual<-list() # keep track of which effective contacts proposed by the specific individual lead to generations
  competition<-list() #keep track of all the effective contact intervals potential infectors propose to a specific susceptible
  for (i in 1:n){
    effective.contacts[[i]]<-0
    generation.individual[[i]]<-NA
    competition[[i]]<-0
  }
  col.name<-c("infected","time.of.infection","infector")
  dimnames(status.matrix)<-list(NULL,col.name)
  status.matrix[,1]<-0 #all the population is initially susceptible
  
  recovery.vector<-rep(NA,n) #vector giving the recovery times
  infectives<-rep(0,n)
  current.time<-0
  rcv<-rep(Inf,n)
  
  # first infected randomly chosen in the population
  first<-sample(1:n, 1)
  generation$Infector[1]<-first
  status.matrix[first,1] <- 1
  status.matrix[first,2] <- 0
  recovery.vector[first]<-current.time+mu
  rcv[first]<-recovery.vector[first]
  infectives[first]<-1
  time.events<-matrix(NA,1,3)
  time.events[1,]<-c(current.time,1,first)
  
  contact.time<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n),"pr.infectee"=rep(NA,n))   #matrix containing the proposed time of the next possible infectious contact (first colum) 
  # with a randomly selected individual (second column)
  proposed.individual<-0
  index.contact<-rep(0,n) # index that selects the individuals for whom a proposed contact is drawn.
  index.contact[first]<-1
  int.time<-0
  temp.contact.time<-0
  T_g<-0
  recovered<-0
  GT<-NULL
  
  #When only the first pathogen is present
  while(sum(infectives, na.rm = TRUE) > 0){ #while there are still infectives
    #Phase 1: individuals that has to, propose a new social contact
    for (i in which(index.contact==1)){ # for all the individuals that has to propose a global contact
      temp.contact.time<-rexp(1,lambda)+current.time# I generate the next interarrival time for individual i
      index.contact[i]<-0
      if (temp.contact.time<recovery.vector[i]){
        contact.time$pr.infectee[i] <-sample(setdiff(1:n,i),1)
        contact.time$pr.ctc[i]<-temp.contact.time # If it is an infectious contact I keep track of the individual and of the time
      }
    }
    
    #Phase 2: identify the next event: possible infection, recovery or the start of the new pathogen infection
    ifelse(length(which(is.na(contact.time$pr.ctc)==FALSE))>0,T_g<-min(contact.time$pr.ctc, na.rm = T),T_g<-Inf)
    R_a<-min(recovery.vector, na.rm = T)
    
    # Phase 2.2 - next event is a proposed infection with pathogen 1
    if (T_g<R_a){
      current.time<-T_g
      infector<-which(contact.time$pr.ctc ==T_g)
      infectee<-sample(setdiff(1:n,infector),1)
      acceptance.rate<-acc.function.inf.3(t=current.time-status.matrix[infector,2])
      is.eff<-0 #we keep track of a contact that is effective
      if (runif(1)<acceptance.rate){ #the contact is effective
        effective.contacts[[infector]]<-c(effective.contacts[[infector]],current.time)
        is.eff<-1
      }
      
      #competition check
      if (is.eff==1 & ((is.na(status.matrix[infectee,1]))| (!(is.na(status.matrix[infectee,1])) & (status.matrix[infectee,1]==1)))){ #if the contact is effective but made with a recovered or a currently infectious individual
        #      if ((status.matrix[infectee,2]-status.matrix[infector,2])<mu & (status.matrix[infector,2]<status.matrix[infectee,2])){ #we check whether the selected infector could be capable of infecting the selected infectee, i.e. if it was infectious when the infectee susceptible
        if ((status.matrix[infectee,2]< rcv[infector]) & (status.matrix[infector,2]<status.matrix[infectee,2])){ #we check whether the selected infector could be capable of infecting the selected infectee, i.e. if it was infectious when the infectee susceptible          
          competition[[infectee]]<-c(competition[[infectee]],(current.time-status.matrix[infector,2]))            
        }
      }
      
      if (!(is.na(status.matrix[infectee,1])) & (status.matrix[infectee,1]==0) & (is.eff==1)){ #if the contact is effective and the contacted individual is susceptible there is transmission
        status.matrix[infectee,1]<-1
        infectives[infectee]<-1
        status.matrix[infectee,2]<-current.time
        status.matrix[infectee,3]<-infector
        recovery.vector[infectee]<-current.time+mu
        rcv[infectee]<-recovery.vector[infectee]
        index.contact[infectee]<-1
        index.contact[infector]<-1
        GT<-c(GT,(current.time-status.matrix[infector,2]))
        effective.contacts[[infectee]][1]<-current.time #the first element of the effective contact list for a specific individual is his/her infection time
        time.events<-rbind(time.events,c(current.time,1,infectee))
        competition[[infectee]]<-c(current.time, (current.time-status.matrix[infector,2]) )# the first element of the competition list is the generation time that infect the individual
        generation.individual[[infector]]<-c(generation.individual[[infector]],(which(effective.contacts[[infector]]==current.time)-1)) #the effective contact that lead to generation (-1 because the first element of EK is the time of infection)
        generation<-rbind(generation,c(current.time,infector,(current.time-status.matrix[infector,2]),(which(effective.contacts[[infector]]==current.time)-1)))
      }else{
        index.contact[infector]<-1
      }
      contact.time$pr.ctc[infector]<-NA
      #Phase 2.3 a recovery occurs  
    }else{
      current.time<-R_a
      recovered<-which(recovery.vector==R_a)
      recovery.vector[recovered]<-NA
      status.matrix[recovered,1]<-NA
      time.events<-rbind(time.events,c(current.time,0.1,recovered))
      contact.time[recovered,2:5]<-rep(NA,4)
      infectives[recovered]<-NA
      effective.contacts[[recovered]]<-c(effective.contacts[[recovered]],current.time) # the last value of effective contact is the day of recover
    }
  }
  timev.name<-c("time","event","who")
  dimnames(time.events)<-list(NULL,timev.name)
  
  return(list(time.events=time.events, status.matrix=status.matrix, generation=generation, effective.contacts=effective.contacts, generation.individual=generation.individual, competition=competition, GT=GT))
}

sim.comp.depl.hom.bc<-function(mu,lambda,n){
  lambda<-lambda*5
  status.matrix <- matrix(NA,nrow = n,ncol = 3) 
  generation<-data.frame("CurrentTime"=0, "Infector"=0, "GenerationLength"=0,"EffectiveContact"=0)
  effective.contacts<-list() #keep track of the effective contacts each infectious individual makes
  generation.individual<-list() # keep track of which effective contacts proposed by the specific individual lead to generations
  competition<-list() #keep track of all the effective contact intervals potential infectors propose to a specific susceptible
  for (i in 1:n){
    effective.contacts[[i]]<-0
    generation.individual[[i]]<-NA
    competition[[i]]<-0
  }
  col.name<-c("infected","time.of.infection","infector")
  dimnames(status.matrix)<-list(NULL,col.name)
  status.matrix[,1]<-0 #all the population is initially susceptible
  
  recovery.vector<-rep(NA,n) #vector giving the recovery times
  infectives<-rep(0,n)
  current.time<-0
  rcv<-rep(Inf,n)
  
  # first infected randomly chosen in the population
  first<-sample(1:n, 1)
  generation$Infector[1]<-first
  status.matrix[first,1] <- 1
  status.matrix[first,2] <- 0
  recovery.vector[first]<-current.time+mu
  rcv[first]<-recovery.vector[first]
  infectives[first]<-1
  time.events<-matrix(NA,1,3)
  time.events[1,]<-c(current.time,1,first)
  
  contact.time<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n),"pr.infectee"=rep(NA,n))   #matrix containing the proposed time of the next possible infectious contact (first colum) 
  # with a randomly selected individual (second column)
  proposed.individual<-0
  index.contact<-rep(0,n) # index that selects the individuals for whom a proposed contact is drawn.
  index.contact[first]<-1
  int.time<-0
  temp.contact.time<-0
  T_g<-0
  recovered<-0
  GT<-NULL
  
  #When only the first pathogen is present
  while(sum(infectives, na.rm = TRUE) > 0){ #while there are still infectives
    #Phase 1: individuals that has to, propose a new social contact
    for (i in which(index.contact==1)){ # for all the individuals that has to propose a global contact
      temp.contact.time<-rexp(1,lambda)+current.time# I generate the next interarrival time for individual i
      index.contact[i]<-0
      if (temp.contact.time<recovery.vector[i]){
        contact.time$pr.infectee[i] <-sample(setdiff(1:n,i),1)
        contact.time$pr.ctc[i]<-temp.contact.time # If it is an infectious contact I keep track of the individual and of the time
      }
    }
    
    #Phase 2: identify the next event: possible infection, recovery or the start of the new pathogen infection
    ifelse(length(which(is.na(contact.time$pr.ctc)==FALSE))>0,T_g<-min(contact.time$pr.ctc, na.rm = T),T_g<-Inf)
    R_a<-min(recovery.vector, na.rm = T)
    
    # Phase 2.2 - next event is a proposed infection with pathogen 1
    if (T_g<R_a){
      current.time<-T_g
      infector<-which(contact.time$pr.ctc ==T_g)
      infectee<-sample(setdiff(1:n,infector),1)
      acceptance.rate<-acc.function.bc(t=current.time-status.matrix[infector,2])
      is.eff<-0 #we keep track of a contact that is effective
      if (runif(1)<acceptance.rate){ #the contact is effective
        effective.contacts[[infector]]<-c(effective.contacts[[infector]],current.time)
        is.eff<-1
      }
      
      #competition check
      if (is.eff==1 & ((is.na(status.matrix[infectee,1]))| (!(is.na(status.matrix[infectee,1])) & (status.matrix[infectee,1]==1)))){ #if the contact is effective but made with a recovered or a currently infectious individual
        #      if ((status.matrix[infectee,2]-status.matrix[infector,2])<mu & (status.matrix[infector,2]<status.matrix[infectee,2])){ #we check whether the selected infector could be capable of infecting the selected infectee, i.e. if it was infectious when the infectee susceptible
        if ((status.matrix[infectee,2]< rcv[infector]) & (status.matrix[infector,2]<status.matrix[infectee,2])){ #we check whether the selected infector could be capable of infecting the selected infectee, i.e. if it was infectious when the infectee susceptible          
          competition[[infectee]]<-c(competition[[infectee]],(current.time-status.matrix[infector,2]))            
        }
      }
      
      if (!(is.na(status.matrix[infectee,1])) & (status.matrix[infectee,1]==0) & (is.eff==1)){ #if the contact is effective and the contacted individual is susceptible there is transmission
        status.matrix[infectee,1]<-1
        infectives[infectee]<-1
        status.matrix[infectee,2]<-current.time
        status.matrix[infectee,3]<-infector
        recovery.vector[infectee]<-current.time+mu
        rcv[infectee]<-recovery.vector[infectee]
        index.contact[infectee]<-1
        index.contact[infector]<-1
        GT<-c(GT,(current.time-status.matrix[infector,2]))
        effective.contacts[[infectee]][1]<-current.time #the first element of the effective contact list for a specific individual is his/her infection time
        time.events<-rbind(time.events,c(current.time,1,infectee))
        competition[[infectee]]<-c(current.time, (current.time-status.matrix[infector,2]) )# the first element of the competition list is the generation time that infect the individual
        generation.individual[[infector]]<-c(generation.individual[[infector]],(which(effective.contacts[[infector]]==current.time)-1)) #the effective contact that lead to generation (-1 because the first element of EK is the time of infection)
        generation<-rbind(generation,c(current.time,infector,(current.time-status.matrix[infector,2]),(which(effective.contacts[[infector]]==current.time)-1)))
      }else{
        index.contact[infector]<-1
      }
      contact.time$pr.ctc[infector]<-NA
      #Phase 2.3 a recovery occurs  
    }else{
      current.time<-R_a
      recovered<-which(recovery.vector==R_a)
      recovery.vector[recovered]<-NA
      status.matrix[recovered,1]<-NA
      time.events<-rbind(time.events,c(current.time,0.1,recovered))
      contact.time[recovered,2:5]<-rep(NA,4)
      infectives[recovered]<-NA
      effective.contacts[[recovered]]<-c(effective.contacts[[recovered]],current.time) # the last value of effective contact is the day of recover
    }
  }
  timev.name<-c("time","event","who")
  dimnames(time.events)<-list(NULL,timev.name)
  
  return(list(time.events=time.events, status.matrix=status.matrix, generation=generation, effective.contacts=effective.contacts, generation.individual=generation.individual, competition=competition, GT=GT))
}



generation.interval.parallel.fixedpropsuscept.hhscen <- function(R, n, mu, n.InInf, n.Rec){
  
  status.matrix <- matrix(NA,nrow = n,ncol = 3) #matrix containing information about the state of the individuals
  col.name<-c("infected","time.of.infection","infector")
  dimnames(status.matrix)<-list(NULL,col.name)
  status.matrix[,1]<-0 #all the population is initially susceptible
  recovery.vector<-rep(NA,n) #vector of recovery times
  infectives<-rep(0,n) # vector that indicates who is infectious at the current time: 1 infectious 0 non infectious
  current.time<-0
  index.contact<-rep(0,n) # vector that selects the individuals that have to propose a new social contact(global) - 1 yes 0 no
  transmission.parameters<-data.frame("id"=1:n,"q"=rep(NA,n),"total_infectionPeriod"=rep(NA,n), "contact_rate"=R/mu)   #matrix containing the proposed time of the next contact (first colum) and the contact individual (second column)
  time.events<-matrix(NA,1,3)
  lambda<-R/mu
  
  
  #
  #Index cases
  first<-sample(which(status.matrix[,1]==0), n.InInf) #initial case
  for (u in first){
    status.matrix[u,1] <- 1 
    status.matrix[u,2] <- 0
    infectives[u]<-1
    recovery.vector[u]<-current.time+mu # the total length since infection (Exposed+IP) 
    #recovery.vector[first]<-current.time+rexp(1/mu) # the total length since infection (Exposed+IP)
    transmission.parameters$q[u]<-R/lambda
    index.contact[u]<-1 #when 1 individual proposes a contact 
    time.events<-rbind(time.events,c(current.time,1,u))
    transmission.parameters$total_infectionPeriod[u]<-recovery.vector[u]-current.time
  }
  
  
  #Immune
  if (n.Rec>=1){
    status.matrix[sample(which(status.matrix[,1]==0),n.Rec),1]<- -1
  }
  contact.time<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n))   #matrix containing the proposed time of the next contact (first colum) and the contact individual (second column)
  proposed.individual<-0
  int.time<-0
  temp.contact.time<-0
  indiv.prop.ctc<-0
  T_NextCtc<-0
  GT<-NULL
  recovered<-0
  sim.indicator<-1
  new.susceptible<-NULL
  susc.over.time<-data.frame("Time"=current.time, "Susceptibles.prop"=  (length(which(status.matrix[,1]==0))/n) )
  
  fwd.gen.time<-data.frame("InfectionTime"=0, "MeanFwdGenTime"=0)
  ind.fwd.gen.time<-list()
  for (i in 1:n) {
    ind.fwd.gen.time[[i]]<--1    
  }
  ct<-0
  
  #When only the first pathogen is present
  while(sum(infectives, na.rm = TRUE) > 0){ #while there are still infectives
    if (ct==1){break}
    #Phase 1: individuals who has to, propose a new social contact
    for (i in which(index.contact==1) ){ # for all the individuals that has to propose a global contact
      temp.contact.time<-rexp(1,transmission.parameters$contact_rate[i])+current.time# I generate the next interarrival time for individual i
      index.contact[i]<-0
      contact.time$pr.ctc[i]<-temp.contact.time # If it is an infectious contact I keep track of the individual and of the time
    }
    
    #Phase 2: identify the next event: possible infection, someone is quarantined or a recovery
    ifelse(length(which(is.na(contact.time$pr.ctc)==FALSE))>0,T_NextCtc<-min(contact.time$pr.ctc, na.rm = T),T_NextCtc<-Inf) # among all the proposed social contact between houeholds we select the minimum
    R_a<-min(recovery.vector, na.rm = T) # minimum among the recovery times
    
    if (T_NextCtc<R_a){
      NextEvent <-"Contact"
    }else{
      NextEvent <-"Recover"
    }
    
    
    if (NextEvent=="Contact"){
      current.time<-T_NextCtc
      infector<-which(contact.time$pr.ctc ==T_NextCtc) 
      infectee<-sample(setdiff(1:n,infector),1)
      if (status.matrix[infectee,1]==0 ){
        status.matrix[infectee,1]<-1
        status.matrix[infectee,2]<-current.time
        status.matrix[infectee,3]<-infector
        GT<-c(GT,(current.time-status.matrix[infector,2]))
        infectives[infectee]<-1
        recovery.vector[infectee]<-current.time+mu # the total length since infection (Exposed+IP) 
        transmission.parameters$q[infectee]<-R/lambda
        index.contact[infectee]<-1 #when 1 individual proposes a contact 
        time.events<-rbind(time.events,c(current.time,1,infectee))
        transmission.parameters$total_infectionPeriod[infectee]<-recovery.vector[infectee]-current.time
        #save the fwd generation time for the infector
        ind.fwd.gen.time[[infector]]<-c(ind.fwd.gen.time[[infector]], (current.time- status.matrix[infector,2]) )
        #when an individual is infected, a recovered is turned into susceptible, to keep the proportion S/N constant
        ct<-1
      }
      index.contact[infector]<-1
      contact.time$pr.ctc[infector]<-NA
      #Phase 2.3 a recovery occurs
    }
    
    if (NextEvent=="Recover"){
      current.time<-R_a
      for (z in which(recovery.vector==R_a)) {
        recovery.vector[z]<-NA
        status.matrix[z,1]<--1
        contact.time[z,2:3]<-rep(NA,4)
        infectives[z]<-NA
        time.events<-rbind(time.events,c(current.time,0.1,z))
        if (length(ind.fwd.gen.time[[z]])>1){
          fwd.gen.time<-rbind(fwd.gen.time, c(status.matrix[z,2] , mean( ind.fwd.gen.time[[z]] [-1]  ) ))
          ind.fwd.gen.time[[z]]<--1
        }
        
      }
    }
    susc.over.time<-rbind(susc.over.time, c(current.time, (length(which(status.matrix[,1]==0))/n)))
  }
  #When also the other pathogen is present.
  timev.name<-c("time","event","who")
  dimnames(time.events)<-list(NULL,timev.name)
  fwd.gen.time<-fwd.gen.time[-1,]
  
  return(list(time.events=time.events, status.matrix=status.matrix, sim.indicator=sim.indicator, fwd.gen.time=fwd.gen.time, GT=GT, susc.over.time=susc.over.time))
  
}


generation.interval.parallel.vacc <- function(R, n, mu, tv,p){
  
  status.matrix <- matrix(NA,nrow = n,ncol = 3) 
  col.name<-c("infected","time.of.infection","infector")
  dimnames(status.matrix)<-list(NULL,col.name)
  status.matrix[,1]<-0
  nSeed<-1
  
  recovery.vector<-rep(NA,n) #vector giving the recovery times
  
  
  exposed.time<-0
  
  # first infected randomly chosen in the population
  first<-sample(1:n, nSeed)
  status.matrix[first,1] <- 1
  status.matrix[first,2]<-0
  
  
  #recovery.vector[first]<-rexp(nSeed,1/mu)+exposed.time
  recovery.vector[first]<-mu+exposed.time
  
  lambda<-R/mu #hazard of the ICID-Exp
  
  
  current.time<-0
  time.events<-matrix(NA,1,3)
  time.events[1,]<-c(0,1,first[1])
  if (nSeed>1){
    for(i in 2:nSeed){
      time.events<-rbind(time.events,c(0,1,first[i]))       
    }  
  }
  
  proposed.contact.time<-matrix(NA,n,2)   #matrix containing the proposed time of the next possible infectious contact (first colum) 
  # with a randomly selected individual (second column)
  proposed.individual<-0
  index.contact<-rep(0,n) # index that selects the individuals for whom a proposed contact is drawn.
  index.contact[first]<-1
  jump<-0
  int.time<-0
  prob.susc<-matrix(NA,1,2)
  competition<-matrix(NA,1,6)
  competition.susc<-matrix(NA,1,6)
  
  check.var<-rep(0,n)
  cont<-0
  
  #first
  while(sum(status.matrix[,1], na.rm = TRUE) > 0){ #while there are still infectives
    contact.times<-rep(NA,n)
    
    if (current.time>tv & cont==0){ # 1 vaccination shot right after the vaccination time tv
      ty<-round(p*length(which(status.matrix[,1]==0))) # number of individuals are going to be vaccinated
      for (q in 1:ty){
        ll<-sample(which(status.matrix[,1]==0),1,replace = FALSE)
        status.matrix[ll,1]<-NA
      }
      cont<-1
    }
    
    for (i in which(index.contact==1) ){ # for all the individuals that has to propose a new contact
      contact.times[i]<-rexp(1,lambda)+current.time# I generate the next interarrival time for individual i
      if (contact.times[i]<recovery.vector[i]){
        #ki<-proposed.contact.time[which(!is.na(proposed.contact.time[,2])),2]
        #proposed.individual<-sample(setdiff(1:n,c(i,ki)),1)
        proposed.individual<-sample(setdiff(1:n,i),1)
        proposed.contact.time[i,]<-c(contact.times[i],proposed.individual) # If it is an infectious contact I keep track of the individual and of the time
      }
    }
    #Competition: the matrix proposed.contact time contains in the first colum the proposed infectious contact time
    # and in the second column the individual the individual to who the infectious contact is proposed.
    # In this matrix are stored both the proposed contacts made in this time step and the proposed contact 
    # made in other time step by individuals that are still infectious
    
    temp.var<-proposed.contact.time[,2]
    temp.var[is.na(temp.var)]<-0
    
    if (!isTRUE(all.equal(temp.var,check.var))){ #here I check if someone new make a contact; in that case I compute the competition. Can be that noone proposed contact and without checking this the program will run always this part
      vector.proposed.individual<-proposed.contact.time[,2] # vector containing for each individual whether he or she proposes a contact (NA or not) indicating the individual to who this infectious contact is directed
      eff.proposed.individual<-vector.proposed.individual[which(!is.na(vector.proposed.individual))]    
      duplicate.pr.ind<-as.vector(unique(eff.proposed.individual[ which(duplicated(eff.proposed.individual)==TRUE)])) #if there more repetitions means that more than one infectious contact the same individual. In this vector we store individuals affect by multiple source of infections
      multiple.infectors<-0  
      if (length(duplicate.pr.ind)>0){
        for(o in 1:length(duplicate.pr.ind)){
          multiple.infectors<-which(vector.proposed.individual==duplicate.pr.ind[o]) #the index identifying the infectors that are competing
          proposed.generation<-0 
          for(g in 1:length(multiple.infectors)){
            proposed.generation[g]<-proposed.contact.time[multiple.infectors[g],1]-status.matrix[multiple.infectors[g],2] #generations proposed by these infectors
          }
          competition<-rbind(competition, c(duplicate.pr.ind[o],min(proposed.generation),max(proposed.generation),length(multiple.infectors), current.time,status.matrix[duplicate.pr.ind[o],1]))
        }
      }
      check.var<-temp.var
    }              
    
    
    if (length(which(is.na(proposed.contact.time[,1])==FALSE))>0){ #If there is at least one proposed contact
      temp2<-min(proposed.contact.time[,1], na.rm = T)
      index.contact<-rep(0,n)
      current.time<-temp2
      infector<-which(proposed.contact.time[,1]==temp2)
      prob.susc<-rbind(prob.susc,c(current.time,length(which(status.matrix[,1]==0))/n))
      if (length(infector)>1){
        temp.inf<-infector
        for (u in 1:length(temp.inf)){
          infector<-temp.inf[u]
          if (!(is.na(status.matrix[proposed.contact.time[infector,2],1])) & status.matrix[proposed.contact.time[infector,2],1]==0){
            status.matrix[proposed.contact.time[infector,2],1]<-1
            index.contact[proposed.contact.time[infector,2]]<-1
            index.contact[infector]<-1
            status.matrix[proposed.contact.time[infector,2],2]<-current.time
            recovery.vector[proposed.contact.time[infector,2]]<-current.time+mu
            #                        recovery.vector[proposed.contact.time[infector,2]]<-current.time+rexp(1,1/mu)
            status.matrix[proposed.contact.time[infector,2],3]<-infector
            time.events<-rbind(time.events,c(current.time,1,proposed.contact.time[infector,2]))
          }else{
            jump<-jump+1
            index.contact[infector]<-1
          }
        }
      }else{
        if (!(is.na(status.matrix[proposed.contact.time[infector,2],1])) & status.matrix[proposed.contact.time[infector,2],1]==0){
          status.matrix[proposed.contact.time[infector,2],1]<-1
          index.contact[proposed.contact.time[infector,2]]<-1
          index.contact[infector]<-1
          status.matrix[proposed.contact.time[infector,2],2]<-current.time
          recovery.vector[proposed.contact.time[infector,2]]<-current.time+mu
          #recovery.vector[proposed.contact.time[infector,2]]<-current.time+rexp(1,1/mu)
          status.matrix[proposed.contact.time[infector,2],3]<-infector
          time.events<-rbind(time.events,c(current.time,1,proposed.contact.time[infector,2]))
        }else{
          jump<-jump+1
          index.contact[infector]<-1
        }
      }
      proposed.contact.time[which(proposed.contact.time[,1]==temp2),]<-rep(NA,2)
    }else{
      current.time<-min(contact.times, na.rm = T)
    }
    
    for (z in which(status.matrix[,1]==1) ){ #recover individual
      if (current.time>recovery.vector[z]){
        status.matrix[z,1]<-NA
        contact.times[z]<-NA
        time.events<-rbind(time.events,c(recovery.vector[z],0,z))
      }
    }
    #proposed.contact.time[,2]
    
  }
  #
  temp.time.events<-time.events #reorder the time events in a chronological order
  
  if (nSeed>1){
    temp.time.events<-temp.time.events[-(1:nSeed),]
    for (i in 1:(length(temp.time.events[,1])-1)){
      temp<-which(min(temp.time.events[,1])==temp.time.events[,1])
      for (j in 1:length(temp)){
        time.events[nSeed+i+(j-1),]<-temp.time.events[temp[j],] 
        temp.time.events<-temp.time.events[-temp[j],]          
      }
      
    }
    time.events[length(time.events[,1]),]<-temp.time.events
    
  }else{
    temp.time.events<-time.events #reorder the time events in a chronological order
    for (i in 1:(length(time.events[,1])-1)){
      temp<-which(min(temp.time.events[,1])==temp.time.events[,1])
      time.events[i,]<-temp.time.events[temp,]
      temp.time.events<-temp.time.events[-temp,]
    }
    time.events[length(time.events[,1]),]<-temp.time.events
  }
  
  timev.name<-c("time","event","who")
  dimnames(time.events)<-list(NULL,timev.name)
  
  competition.susc<-competition[which(competition[,6]==0),]
  if (is.matrix(competition.susc) && length(competition.susc[,1])>2){ # clean the competition matrix. Delete the same reported value
    temp.compe<-matrix(NA,1,6)
    competition2<-competition.susc
    #competition2<-as.matrix(competition.susc[-1,])
    competition2[is.na(competition.susc[,6]),6]<-2
    temp.compe[1,]<-competition2[1,]
    
    
    for (i in 2:length(competition2[,1])){
      count<-1
      for(j in (i-1):1){
        if((competition2[i,1]==competition2[j,1]) & (competition2[i,2]==competition2[j,2]) & (competition2[i,3]==competition2[j,3]) & (competition2[i,4]==competition2[j,4])){
          count<-0
        }
      }
      if(count==1){
        temp.compe<-rbind(temp.compe,competition2[i,])
      }
    }
    competition.susc<-temp.compe
  }
  comp.name<-c("contacted", "Min","Max", "N-comp.","current.time","status.contacted")
  dimnames(competition)<-list(NULL,comp.name)
  
  
  temp<-0
  del.line<-0
  if(is.matrix(competition.susc)){
    if (length(competition.susc[,1])!=length(unique(competition.susc[,1]))){
      uniq.sucept<-duplicated(competition.susc[,1])
      temp<-competition.susc
      for (b in which(uniq.sucept==TRUE)){
        temp2<-which(competition.susc[,1]==competition.susc[b,1])
        temp3<-which(max(competition.susc[temp2,4])==competition.susc[temp2,4])
        del.line<-c(del.line,setdiff(temp2,temp2[temp3]))
      }
      competition.susc<-competition.susc[-del.line[1:length(del.line)],]
    }
  }
  
  
  # List of output variables
  return(list(time.events=time.events, status.matrix=status.matrix, jump=jump, prob.susc=prob.susc, competition=competition, competition.susc=competition.susc) )
}


sim.comp.depl.vacc<-function(mu,lambda,n, tv,p){
  
  status.matrix <- matrix(NA,nrow = n,ncol = 3) 
  generation<-data.frame("CurrentTime"=0, "Infector"=0, "GenerationLength"=0,"EffectiveContact"=0)
  effective.contacts<-list() #keep track of the effective contacts each infectious individual makes
  generation.individual<-list() # keep track of which effective contacts proposed by the specific individual lead to generations
  competition<-list() #keep track of all the effective contact intervals potential infectors propose to a specific susceptible
  for (i in 1:n){
    effective.contacts[[i]]<-0
    generation.individual[[i]]<-NA
    competition[[i]]<-0
  }
  col.name<-c("infected","time.of.infection","infector")
  dimnames(status.matrix)<-list(NULL,col.name)
  status.matrix[,1]<-0 #all the population is initially susceptible
  
  recovery.vector<-rep(NA,n) #vector giving the recovery times
  infectives<-rep(0,n)
  current.time<-0
  rcv<-rep(Inf,n)
  
  # first infected randomly chosen in the population
  first<-sample(1:n, 1)
  generation$Infector[1]<-first
  status.matrix[first,1] <- 1
  status.matrix[first,2] <- 0
  recovery.vector[first]<-current.time+mu
  rcv[first]<-recovery.vector[first]
  infectives[first]<-1
  time.events<-matrix(NA,1,3)
  time.events[1,]<-c(current.time,1,first)
  
  contact.time<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n),"pr.infectee"=rep(NA,n))   #matrix containing the proposed time of the next possible infectious contact (first colum) 
  # with a randomly selected individual (second column)
  proposed.individual<-0
  index.contact<-rep(0,n) # index that selects the individuals for whom a proposed contact is drawn.
  index.contact[first]<-1
  int.time<-0
  temp.contact.time<-0
  T_g<-0
  recovered<-0
  GT<-NULL
  cont<-0
  
  #When only the first pathogen is present
  while(sum(infectives, na.rm = TRUE) > 0){ #while there are still infectives
    #Phase 1: individuals that has to, propose a new social contact
    for (i in which(index.contact==1)){ # for all the individuals that has to propose a global contact
      temp.contact.time<-rexp(1,lambda)+current.time# I generate the next interarrival time for individual i
      index.contact[i]<-0
      if (temp.contact.time<recovery.vector[i]){
        contact.time$pr.infectee[i] <-sample(setdiff(1:n,i),1)
        contact.time$pr.ctc[i]<-temp.contact.time # If it is an infectious contact I keep track of the individual and of the time
      }
    }
    
    #Phase 2: identify the next event: possible infection, recovery or the start of the new pathogen infection
    ifelse(length(which(is.na(contact.time$pr.ctc)==FALSE))>0,T_g<-min(contact.time$pr.ctc, na.rm = T),T_g<-Inf)
    R_a<-min(recovery.vector, na.rm = T)
    
    
    if (current.time>tv & cont==0){ # 1 vaccination shot right after the vaccination time tv
      ty<-round(p*length(which(status.matrix[,1]==0))) # number of individuals are going to be vaccinated
      for (q in 1:ty){
        ll<-sample(which(status.matrix[,1]==0),1,replace = FALSE)
        status.matrix[ll,1]<--2
        status.matrix[ll,2]<-current.time
      }
      cont<-1
    }
    
    
    # Phase 2.2 - next event is a proposed infection with pathogen 1
    if (T_g<R_a){
      current.time<-T_g
      infector<-which(contact.time$pr.ctc ==T_g)
      infectee<-sample(setdiff(1:n,infector),1)
      acceptance.rate<-acc.function.k(t=current.time-status.matrix[infector,2])
      is.eff<-0 #we keep track of a contact that is effective
      if (runif(1)<acceptance.rate){ #the contact is effective
        effective.contacts[[infector]]<-c(effective.contacts[[infector]],current.time)
        is.eff<-1
      }
      
      #competition check
      if (is.eff==1 & ((is.na(status.matrix[infectee,1]))| (!(is.na(status.matrix[infectee,1])) & (status.matrix[infectee,1]==1)))){ #if the contact is effective but made with a recovered or a currently infectious individual
        #      if ((status.matrix[infectee,2]-status.matrix[infector,2])<mu & (status.matrix[infector,2]<status.matrix[infectee,2])){ #we check whether the selected infector could be capable of infecting the selected infectee, i.e. if it was infectious when the infectee susceptible
        if ((status.matrix[infectee,2]< rcv[infector]) & (status.matrix[infector,2]<status.matrix[infectee,2])){ #we check whether the selected infector could be capable of infecting the selected infectee, i.e. if it was infectious when the infectee susceptible          
          competition[[infectee]]<-c(competition[[infectee]],(current.time-status.matrix[infector,2]))
          if (status.matrix[infectee,1]==-2){competition[[infectee]]<-0}
        }
      }
      
      if (!(is.na(status.matrix[infectee,1])) & (status.matrix[infectee,1]==0) & (is.eff==1)){ #if the contact is effective and the contacted individual is susceptible there is transmission
        status.matrix[infectee,1]<-1
        infectives[infectee]<-1
        status.matrix[infectee,2]<-current.time
        status.matrix[infectee,3]<-infector
        recovery.vector[infectee]<-current.time+mu
        rcv[infectee]<-recovery.vector[infectee]
        index.contact[infectee]<-1
        index.contact[infector]<-1
        GT<-c(GT,(current.time-status.matrix[infector,2]))
        effective.contacts[[infectee]][1]<-current.time #the first element of the effective contact list for a specific individual is his/her infection time
        time.events<-rbind(time.events,c(current.time,1,infectee))
        competition[[infectee]]<-c(current.time, (current.time-status.matrix[infector,2]) )# the first element of the competition list is the generation time that infect the individual
        generation.individual[[infector]]<-c(generation.individual[[infector]],(which(effective.contacts[[infector]]==current.time)-1)) #the effective contact that lead to generation (-1 because the first element of EK is the time of infection)
        generation<-rbind(generation,c(current.time,infector,(current.time-status.matrix[infector,2]),(which(effective.contacts[[infector]]==current.time)-1)))
      }else{
        index.contact[infector]<-1
      }
      contact.time$pr.ctc[infector]<-NA
      #Phase 2.3 a recovery occurs  
    }else{
      current.time<-R_a
      recovered<-which(recovery.vector==R_a)
      recovery.vector[recovered]<-NA
      status.matrix[recovered,1]<-NA
      time.events<-rbind(time.events,c(current.time,0.1,recovered))
      contact.time[recovered,2:5]<-rep(NA,4)
      infectives[recovered]<-NA
      effective.contacts[[recovered]]<-c(effective.contacts[[recovered]],current.time) # the last value of effective contact is the day of recover
    }
  }
  timev.name<-c("time","event","who")
  dimnames(time.events)<-list(NULL,timev.name)
  
  return(list(time.events=time.events, status.matrix=status.matrix, generation=generation, effective.contacts=effective.contacts, generation.individual=generation.individual, competition=competition, GT=GT))
}




sim.comp.depl.2lvlmixing<-function(HH.network, lambda.g, nSeeds,inf.path.h,inf.path.g){
  
  n<-network.size(HH.network)
  hh.id<- HH.network %v% "hh_id"
  status.matrix <- matrix(NA,nrow = n,ncol = 7) #matrix containing information about the state of the individuals
  col.name<-c("infected","time.of.infection","infector", "severity","Exposed", "TimeSymptomOnset","RecoveryTime")
  dimnames(status.matrix)<-list(NULL,col.name)

  
  generation<-data.frame("CurrentTime"=0, "Infector"=0, "GenerationLength"=0,"EffectiveContact"=0)
  effective.contacts<-list() #keep track of the effective contacts each infectious individual makes
  generation.individual<-list() # keep track of which effective contacts proposed by the specific individual lead to generations
  competition.g<-list() #keep track of all the effective contact intervals potential infectors propose to a specific susceptible
  competition.l<-list()
  for (i in 1:n){
    effective.contacts[[i]]<-0
    generation.individual[[i]]<-NA
    competition.g[[i]]<-0
    competition.l[[i]]<-0
  }
  GT<-NULL
  
  status.matrix[,1]<-0 #all the population is initially susceptible
  status.matrix[,6]<-Inf 
  recovery.vector<-rep(Inf,n) #vector giving the recovery times

  events<-data.frame("NextCtc"=Inf, "Recovery"=Inf)

  infectives<-rep(0,n) # vector that indicates who is infectious at the current time: 1 infectious 0 non infectious
  current.time<-0
  index.contact.within<-rep(0,n) # vector that selects the individuals that have to propose a new social contact(global) - 1 yes 0 no
  index.contact.between<-rep(0,n) # vector that selects the individuals that have to propose a new social contact(global) - 1 yes 0 no
  time.events<-matrix(NA,1,3)
  rcv<-rep(Inf,n)
  
  #transmission parameter dataframe: each line is an individual, the first colum is the transmsission coeffficient and the third the length of IP (needed to re-scale Viral load curve)
  transmission.parameters<-data.frame("id"=1:n,"qh"=rep(NA,n),"qg"=rep(NA,n),"contact_rate_within"=rep(NA,n),"contact_rate_between"=lambda.g, "susceptibility"=rep(1,n))   #matrix containing the proposed time of the next contact (first colum) and the contact individual (second column)
  
  for (j in 1:n){
    transmission.parameters$contact_rate_within[j]<-length(get.neighborhood(HH.network,j))
  }
  
  #Proportion of immune
  if (prop.immune>0){
    status.matrix[sample(1:n,round(prop.immune*n)),1]<--2
  }
  

  contact.time.within<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n),"pr.infectee"=rep(NA,n))   #matrix containing the proposed time of the next contact (first colum) and the contact individual (second column)
  contact.time.between<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n),"pr.infectee"=rep(NA,n))   #matrix containing the proposed time of the next contact (first colum) and the contact individual (second column)
  
  # first infected: randomly chosen in the population (among the susceptibles)
  first.cases<-sample(which(status.matrix.1[,1]==0),nSeeds)
  
  for (j in first.cases){
    first<-j
    status.matrix[first,1] <- 1 
    status.matrix[first,2] <- 0
    status.matrix[first,5] <- current.time+exposed.time()
    status.matrix[first,7]<-status.matrix[first,5]+infectious.period()
    recovery.vector[first]<-status.matrix[first,7]
    transmission.parameters$qh[first]<-inf.path.h #A single q parameter for everyone
    transmission.parameters$qg[first]<-inf.path.g #A single q parameter for everyone
    time.events<-rbind(time.events,c(current.time,1,first))
    infectives[first]<-1
    contact.time.within$pr.ctc[first]<-rexp(1,transmission.parameters$contact_rate_within[first])+current.time+status.matrix.1[first,5]# I generate the next interarrival time for individual i
    contact.time.between$pr.ctc[first]<-rexp(1,transmission.parameters$contact_rate_between[first])+current.time+status.matrix.1[first,5] # I generate the next interarrival time for individual i
  }
  
  proposed.individual<-0
  temp.contact.time<-0
  indiv.prop.ctc<-0
  recovered<-0
  err<-0
  
  #When only the first pathogen is present
  while((sum(infectives))>0 | current.time<t2){ #while there are still infectives
    #Phase 1: individuals that has to, propose a new social contac
    
    for (i in which(index.contact.within==1) ){ # for all the individuals that has to propose a global contact
      contact.time.within$pr.ctc[i]<-rexp(1,transmission.parameters$contact_rate_within[i])+current.time# I generate the next interarrival time for individual i
      index.contact.within[i]<-0
    }
    for (i in which(index.contact.between==1) ){ # for all the individuals that has to propose a global contact
      contact.time.between$pr.ctc[i]<-rexp(1,transmission.parameters$contact_rate_between[i])+current.time# I generate the next interarrival time for individual i
      index.contact.between[i]<-0
    }
    
    contact.time.overall<-c(contact.time.within$pr.ctc, contact.time.between$pr.ctc) #overall contact times    
    #Phase 2: identify the next event: possible infection, recovery or the start of the new pathogen infection
    ifelse(length(which(is.na(contact.time.overall)==FALSE))>0,events$NextCtc<-min(contact.time.overall, na.rm = T),events$NextCtc<-Inf) # among all the proposed social contact between houeholds we select the minimum
    ifelse(length(which(is.na(recovery.vector)==FALSE))>0,events$Recovery<-min(recovery.vector, na.rm = T),events$Recovery<-Inf) # among all the proposed social contact between houeholds we select the minimum
    
    next.evts<-colnames(events)[which(min(events)==events)]
    if (length(next.evts)>1){
      next.evts<-sample(colnames(events)[which(min(events)==events)],1)
    }
    if (next.evts=="NextCtc"){
      current.time<-events$NextCtc
      if (length(min(contact.time.overall, na.rm = T))>1){ #when two contacts happen at the same time
        selected.ctc<-sample(which(contact.time.overall==current.time),1) 
        if (selected.ctc!=n & selected.ctc!=2*n){
          infector<- selected.ctc %% n
        }else{
          infector<- n
        }
        if (selected.ctc<=n){
          infectee.pool<-get.neighborhood(HH.network,infector)
          if (length(infectee.pool)>1){
            infectee<-sample(infectee.pool,1) #pick a random individual in the class
          }
          if(length(infectee.pool)==1){
            infectee<-infectee.pool #pick a random individual in the class
          }
          if(length(infectee.pool)==0){
            infectee<-infector #just a trick to have acceptance rate 0 (infector is not susceptible)
          }
          index.contact.within[infector]<-1
          contact.time.within$pr.ctc[infector]<-NA
          ctc<-"hh"
        }else{
          hh.members.temp<-which(hh.id==hh.id[infector])
          infectee.pool<- setdiff(1:n,hh.members.temp) #individuals not in the same class
          infectee<-sample(infectee.pool,1) #pick a random individual not in the class (and not a teacher)
          index.contact.between[infector]<-1
          contact.time.between$pr.ctc[infector]<-NA
          ctc<-"g"
        }
      }else{
        if (length(which(events$NextCtc==contact.time.within$pr.ctc))>0){ #if it is a within contact
          infector<-which(contact.time.within$pr.ctc ==events$NextCtc)
          infectee.pool<-get.neighborhood(HH.network,infector)
          if (length(infectee.pool)>1){
            infectee<-sample(infectee.pool,1) #pick a random individual in the class
          }
          if(length(infectee.pool)==1){
            infectee<-infectee.pool #pick a random individual in the class
          }
          if(length(infectee.pool)==0){
            infectee<-infector #just a trick to have acceptance rate 0 (infector is not susceptible)
          }
          index.contact.within[infector]<-1
          contact.time.within$pr.ctc[infector]<-NA
          ctc<-"hh"
        }else{
          infector<-which(contact.time.between$pr.ctc == events$NextCtc) 
          hh.members.temp<-which(hh.id==hh.id[infector])
          infectee.pool<- setdiff(1:n,hh.members.temp) #individuals not in the same class
          infectee<-sample(infectee.pool,1) #pick a random individual not in the class (and not a teacher)
          index.contact.between[infector]<-1
          contact.time.between$pr.ctc[infector]<-NA
          ctc<-"g"
        }
      }
      
      # compute the long and short interaction terms for pathogen.1
   
      ifelse(ctc=="g",q<-transmission.parameters$qg[infector],q<-transmission.parameters$qh[infector])
      acc.rate<-InfMeasure(t=current.time-status.matrix[infector,5])*q
      if (acc.rate.1>1){err<-err+1}
      is.effective<-0
      if (runif(1)<acc.rate){
        effective.contacts[[infector]]<-c(effective.contacts[[infector]],current.time)
        is.effective<-1
      }
      
      #Competition check
      if ( (is.eff==1 & (status.matrix[infectee,1]==-1)) |  (is.eff==1 & (status.matrix[infectee,1]==1))){ #if the contact is effective but made with a recovered or a currently infectious individual
        if ((status.matrix[infectee,2]< status.matrix[infector,7]) & (status.matrix[infector,2]<status.matrix[infectee,2])){ #we check whether the selected infector could be capable of infecting the selected infectee, i.e. if it was infectious when the infectee susceptible          
          
          if (ctc=="g"){
            competition.g[[infectee]]<-c(competition.g[[infectee]],(current.time-status.matrix[infector,2]))            
          }else{
            competition.l[[infectee]]<-c(competition.l[[infectee]],(current.time-status.matrix[infector,2]))            
            
          }
          
        }
      }
      
      if (status.matrix[infectee,1]==0 & is.effective==1){
        status.matrix[infectee,1] <- 1 
        status.matrix[infectee,2] <- current.time
        status.matrix[infectee,3] <- infector
        status.matrix[infectee,5] <- current.time+exposed.time.1()
        status.matrix[infectee,7]<-status.matrix[infectee,5]+infectious.period()
        recovery.vector[infectee]<-status.matrix[infectee,7]
        transmission.parameters$qh[infectee]<-inf.path.h #A single q parameter for everyone
        transmission.parameters$qg[infectee]<-inf.path.g #A single q parameter for everyone
        time.events<-rbind(time.events,c(current.time,1,infectee))
        GT<-c(GT,(current.time-status.matrix[infector,2]))
        effective.contacts[[infectee]][1]<-current.time #the first element of the effective contact list for a specific individual is his/her infection time
        competition.g[[infectee]]<-c(current.time, (current.time-status.matrix[infector,2]) )# the first element of the competition list is the generation time that infect the individual
        competition.g[[infectee]]<-c(current.time, (current.time-status.matrix[infector,2]) )# the first element of the competition list is the generation time that infect the individual
        generation.individual[[infector]]<-c(generation.individual[[infector]],(which(effective.contacts[[infector]]==current.time)-1)) #the effective contact that lead to generation (-1 because the first element of EK is the time of infection)
        generation<-rbind(generation,c(current.time,infector,(current.time-status.matrix[infector,2]),(which(effective.contacts[[infector]]==current.time)-1)))
      }
      
      
      # compute the long and short interaction terms for pathogen.2
    }
    
    if (next.evts=="Recovery"){
      current.time<-events$Recovery
      temp.recovered<-which(recovery.vector.overall==events$Recovery)
      for (recovered in temp.recovered){
            status.matrix[recovered,1]<--1
            recovery.vector[recovered]<-Inf
            time.events<-rbind(time.events,c(current.time,-1,recovered))
            infectives[recovered]<-0
            contact.time.between$pr.ctc[recovered]<-NA
            contact.time.within$pr.ctc[recovered]<-NA
            index.contact.within[recovered]<-0
            index.contact.between[recovered]<-0
              }
    }
  }
  
  #When also the other pathogen is present.
  time.events<-time.events[-1,]
  timev.name<-c("time","event","who")
  dimnames(time.events)<-list(NULL,timev.name)
  
  first.cases<-which(status.matrix[,2]==0)
  for (o in first.cases){
    temp.sec.cases<-NULL
    ifelse(length(which(status.matrix[,3]==o)>0),temp.sec.cases<-c(temp.sec.cases,length(which(status.matrix[,3]==o))),temp.sec.cases<-c(temp.sec.cases,0))
  }
  Rt1<-mean(temp.sec.cases)
  
  C1<-nSeeds.1
  Y1<-nSeeds.1
  last.day<-round(max(time.events[,1]))
  
  for (i in 1:last.day){
    temp.time<-setdiff(which(time.events[,1]>i),which(time.events[,1]>i+1))
    temp.inf.1<-c(which(time.events[temp.time,2]==1.0))
    temp.time.1<-setdiff(1:length(time.events[,1]),which(time.events[,1]>i+1))
    C1<- c(C1,length((which(time.events[temp.time.1,2]==1.0)))+length((which(time.events[temp.time.1,2]==1.1)))-length((which(time.events[temp.time.1,2]==-1))))
    Y1<- c(Y1,length(temp.inf.1))

    if (length(temp.inf.1)>0){
      newly.infected<-time.events[temp.time[temp.inf.1],3]
      temp.sec.cases<-NULL
      for (k in newly.infected) {
        ifelse(length(which(status.matrix[,3]==k)>0),temp.sec.cases<-c(temp.sec.cases,length(which(status.matrix[,3]==k))),temp.sec.cases<-c(temp.sec.cases,0))
      }
      Rt1<-c(Rt1,mean(temp.sec.cases))
    }else{
      Rt1<-c(Rt1,NA)
    }
  }
  epi.details<-data.frame("Days"=0:last.day, "Incidence1"=Y1, "Prevalence1"=C1,"Rt1"=Rt1)
  
  FinalSize<-data.frame("FinalSize1"=length(which(status.matrix.1==-1)),"FinalSize2"=length(which(status.matrix.2==-1)))
  return(list(time.events=time.events, status.matrix.1=status.matrix.1,epi.details=epi.details, FinalSize=FinalSize,generation=generation, effective.contacts=effective.contacts, generation.individual=generation.individual, competition.l=competition.l,competition.g=competition.g, GT=GT))
}









