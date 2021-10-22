
SIR_pomp_function<-function(incubationT,delta,r0,reportingRate,L,pop,nSim,paso_tiempo,tiempos,semilla){

  I_t0<-0.0*pop #initial number of infecions in the units
  R_t0<-0.0*pop #initial number of recover people in the units
  S_t0<-pop-I_t0-R_t0 # initial number of susceptibles
  A_beta<-rep(r0/incubationT,L)# paramenter related to transmission rate
  B_beta<-rep(delta,L)## paramenter related to transmission rate

  global_params=c(mu_H = 3.68e-05, mu_EI = 0.5984005,gamma = 1/incubationT, rho = reportingRate,
                  omega = (2*pi/365), phi = 1.159567, sigma_P = 0.28899, sigma_M = 0.1409216)
  params = c(I_t0, R_t0, S_t0, A_beta, B_beta, global_params, L)
  names(params)=paramnames
  
  #### simulation ####

  sim_data = simulate(
    format="data.frame",times=tiempos, t0=0,covar=covariate_table(sparks,times="time"),
    nsim = nSim, seed = semilla,
    rprocess = euler(rproc,delta.t = paso_tiempo),
    params = params, paramnames = paramnames,
    statenames = statenames, obsnames = obsnames, accumvars = c("C"),
    rinit = rInitVic,
    rmeas = rmeas)
  return(sim_data)
}

### L style ###
format_trasf<-function(nSim,data,tf){
  output_fromat<-array(NA,dim=c(tf,nSim))
  auxCases<-data[,"C"]
  id<-as.numeric(data[,".id"])
  for(n in 1:nSim){
    
    who<-which(id==n)
    output_fromat[,n]<-auxCases[who]
    rm(who)
  }
  return(output_fromat)
}


# function to go from daily to monthly #
FromDaysToMoth_function<-function(nYear,CasesDays,nSim,month0){
  Days_year<-c(31,28,31,30,31,30,31,31,30,31,30,31)
  Days_First<-Days_year[c(month0:length(Days_year))]
  aux<-c(Days_First,Days_year)
  aux<-aux[c(1:12)]
  
  DaysMonthC<-rep(aux, nYear)
  if(nSim>1){
    nUnits<-nSim
 
    aux_ObsSimCases<-array(NA,dim=c(length(DaysMonthC),nUnits))
  
    sup<-0
    for(d in 1:length(DaysMonthC)){
      inf<-1+sup
      sup<-sup+DaysMonthC[d]
      aux_ObsSimCases[d,]<-apply(CasesDays[c(inf:sup),],2,sum)

      
    }
  }
  if(nSim==1){
    aux_ObsSimCases<-rep(NA,length(DaysMonthC))
    sup<-0
    for(d in 1:length(DaysMonthC)){
      inf<-1+sup
      sup<-sup+DaysMonthC[d]
      
      aux_ObsSimCases[d]<-sum(CasesDays[c(inf:sup)])
      #s_input[d,p]<-sum(a[c(inf:sup)])
      
    }
  }
  return(aux_ObsSimCases)
}

Select_pop<-function(popG){
  
  aux<-which(is.na(popG)==TRUE)
  output_pop<-popG
  if(length(aux)>0)
    output_pop<-popG[-aux]
  return(output_pop)
}

# uniformly distribution from month to days #
FromMonthToDays_Uniform<-function(inputVector,TimeSaco,Years_FromJan){

	daysPerMonth<-rep(c(31,28,31,30,31,30,31,31,30,31,30,31), length(Years_FromJan))
	output_UniformPerDay<-rep(NA,sum(daysPerMonth))

	i_m<-1
	i_M<-0
	for(m in 1:length(daysPerMonth)) {
		i_M<-i_M+daysPerMonth[m]
		output_UniformPerDay[c(i_m:i_M)]<-rep(inputVector[m]/daysPerMonth[m],daysPerMonth[m])
		i_m<-i_M+1
	}
        if(TimeSaco>1)
		output_UniformPerDay<-output_UniformPerDay[-c(1:TimeSaco)]
	return(output_UniformPerDay)


}

# function to compute the arrival time #
ArrivalTimeFunction<-function(rho,infMatrix,recuMatrix,interval,nSim,month0){
  nYears<-length(interval)-1
  tf<-nrow(infMatrix)
  outputArrT<-array(NA,dim=c(ncol(infMatrix),nYears),
                    dimnames=list(NULL,paste0("Season ",as.character(c(1:nYears)))))
  
  for(j in 1:ncol(infMatrix)){
    #aux1<-mapply(rpois,lambda=infMatrix[,j],MoreArgs = list(n=1))
    deltaInfMatrix<-infMatrix[-1,j]-infMatrix[-tf,j] + recuMatrix[-1,j] - recuMatrix[-tf,j]
    deltaInfMatrix[which(deltaInfMatrix<0)]<-0
    #aux1<-mapply(rbinom,size=infMatrix[,j],MoreArgs = list(n=1,prob=rho))
    aux1<-mapply(rbinom,size=deltaInfMatrix,MoreArgs = list(n=1,prob=rho))
    for(i in 1:nYears){
      aux2<-which(FromDaysToMoth_function(nYears,aux1,nSim,month0)[c((interval[i]+1):interval[(i+1)])]>0)
      if(length(aux2)>0)
        outputArrT[j,i]<-min(aux2)
    }
  }
  return(outputArrT)
}

# TheoTotSpark_function<-function(TotCases,Tau){
#  
#   m<-0.624
#   b<-exp(1.854)/Tau
#   output_f<-((TotCases*Tau)^m)*b
#   output_f/Tau
#   return(output_f)
#   #output_f<-0.624*log(TotCases)+1.854
#   #return(exp(output_f))
# }

# TheoSparkPerUnit_function<-function(TotCases,popU,Tau){
# 
#   a<-5.57e-6#1e-6
#   b<-1.17e-3#8.3e-4
#   c<-0.57#0.89
#   d<-0.767#0.733
#   k<- -1.44e-4 #-1.64e-4
#   outputTheoSpark<-a*exp(b*popU)*(popU^c)*(TotCases*Tau)^(d+k*popU)
#   outputTheoSpark<-outputTheoSpark/Tau
#   return(outputTheoSpark)
# }

# spark_per_group_function<-function(TotCaseDay,PopMatrix,Stot,Tau){
#   output_spark_per_group<-array(NA,dim=c(length(TotCaseDay),ncol(PopMatrix)))
#   for(tt in 1:length(TotCaseDay)){
#     lamU<-apply(PopMatrix,2,TheoSparkPerUnit_function,TotCases=TotCaseDay[tt],Tau)
#     lamG<-apply(lamU,2,sum,na.rm=TRUE)
#     output_spark_per_group[tt,]<-rmultinom(n=1,size=Stot[tt],prob=lamG/sum(lamG))
#     if(TotCasesDay[tt]==0){
#       output_spark_per_group[tt,]<-rep(NA,ncol(output_spark_per_group))
#       print("0 total Cases in the city")
#     }
#   }
#   
#   return(output_spark_per_group)
# }


# probSparkPerGroup_function<-function(TotCasesDay,PopMatrix) {
#   if(TotCasesDay==0){
#     print("0 total Cases in the city")
#     prob_g<-rep(NA,length(TotCases))
#   }
#   else{
#     lamU<-TheoTotSpark_function(TotalCasesDay,PopMatrix)
#     prob_g<-apply(lamU,1,sum,na.rm=TRUE)/TheoTotSpark_function(TotCasesDay)
# 
#   }
#   return(prob_g)
# }



# sparkRate_per_group_function<-function(SparksProb,PopMatrix,Stot){
#   output_spark_per_group<-array(NA,dim=dim(SparksProb))
#   Laux<-NULL
#   for(q in 1:12){
#     Laux<-c(Laux,length(which( is.na(PopMatrix[,q])==FALSE)))
#   }
#   for(tt in 2:(nrow(SparksProb)-1)){
#  
#     output_spark_per_group[tt,]<-rmultinom(n=1,size=Stot[tt],prob=SparksProb[tt,])/Laux
#    
#   }
#   
#   return(output_spark_per_group)
# }



