####################################################################################################################
## loading libraries packages
library("pomp")#, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5") # pomp2
# library("plyr")#, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
# library("dplyr")#, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
# library("ggplot2")#, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
# library("reshape2")#, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
# library("lubridate")#, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
#setwd("/home/victoria/Documents/SpatialSparkDengueModel_SIR/Pomp-Version")

nG_vector<-c(12,25,50,100) # number of partitions (based on population quantiles)
Years_FromJan<-as.character(c(2010:2014))
nSim<-1  #number of simulations within pomp. Keep it value in 1. The repetitions are given by Repite variable
paso_tiempo<-1/12 #12 Time step --> time step of 2 hours
tiempos<-c(0:800) #duration of the simulation (days)
interval<-c(0,12,24) #monthes of season. 
TimeSaco<-638 #638 Original #be careful with this. It must agree with t0 in the sinusoidal function beta (see Cisneppt)
month0<-10 # initial month of the dengue season is October
source("Functions_SIR_SparkModel_Prob.R") ## loading functions

# parameters that are varing
incubationT<-17#c(17,12,7)
delta<-c(0.8,0.6,0.4,0.2)
r0<-1.15#c(1.225,1.15, 1.075, 0)
rho<-c(0.1,0.3,0.5)

#the number of times we repeat the simulations are done outside Pomp for memory reasons
Repite<-10 # number of times we repeat the stochastic simulations

# names of files to storage the outputs #
namesPeaks<-"PeakRatio_incT17_r01.15_delta0.2_nQXX.Rdata"
namesCases<-"CasesMonth_incT17_r01.15_delta0.2_nQXX.Rdata"
namesArrT<-"ArrT_incT17_r01.15_delta0.2_nQXX.Rdata"

for(gg in 1:length(nG_vector)){
  nG<-nG_vector[gg]
  PopMatrix<-read.table(file=gsub("XX",as.character(nG),"PopMatrixNoMaxLimit_nGXX.dat"), sep=" ", header = TRUE) 
  nSparksPerUnit<-read.table(file=gsub("XX",as.character(nG),"nSparksMonth_nGXX.dat"),sep = " ",header = TRUE)
  meanPop<-round(apply(PopMatrix,2,mean, na.rm=TRUE),0)
  
  ### variables to storage things of interest #
  peakRatio<-array(NA,dim=c(Repite, nG,length(delta),length(rho)),
                   dimnames=list(NULL,names(meanPop), 
                                 paste0("delta=",as.character(delta)), paste0("rho=",as.character(rho)) ))
  CasesMonth<-array(NA,dim=c(Repite,nG,24,length(delta),length(rho)),
                    dimnames = list(NULL,names(meanPop),paste0("Moth=",as.character(c(1:24))),
                                    paste0("delta=",as.character(delta)), paste0("rho=",as.character(rho)) ))
  
  ArrT<-array(NA,dim=c(Repite,nrow(PopMatrix),(length(interval)-1),nG,length(delta),length(rho)),
              dimnames=list(paste0("Repet_",as.character(1:Repite)),NULL, 
                            paste0("Season_",as.character(1:2)), paste0("PopGroup_",as.character(1:nG)),
                            paste0("delta=",as.character(delta)), paste0("rho=",as.character(rho)) ))

  semilla<-12346 # seed for the simulaitons
  for(r in 1:Repite){
    print(paste0("r=",as.character(r)))
    for(q in 1:nG){
      print(names(meanPop[q])) # population of units belonging to group q
      DailySparkPerGroup<-FromMonthToDays_Uniform(nSparksPerUnit[,q],TimeSaco,Years_FromJan ) # converts monthly sparks to daily (uniformly. that is, dailySparks = monthlySparks / number_days_in_the_month
      
      sparks<-data.frame(time=c(0:(length(DailySparkPerGroup)-1)),
                         sp_rate=DailySparkPerGroup) # covarite table for pomp object input
      
      pop<-Select_pop(PopMatrix[,q]) # population of the units that belong to group q
      L<-length(pop) # number of units in group q
      source("Csnippet_SIR_SparkSpatialModel.R") # it has the process model, names of parameters and variables.
      
      for(d in 1:length(delta)){ # loop for parameter delta
        for(rho_i in 1:length(rho)){ # loop for reporting rate values
          
          sim_data<-SIR_pomp_function(incubationT,delta[d],r0,rho[rho_i],L,pop,nSim,paso_tiempo,tiempos,semilla) # SIR simulation of all units in q group
          
          ArrT[r,c(1:L),,q,d,rho_i]<-ArrivalTimeFunction(rho[rho_i],infMatrix = sim_data[,paste0("I",as.character(1:L))],
                                                         recuMatrix = sim_data[,paste0("R",as.character(1:L))],interval,
                                                         nSim,month0) # computation of the arrival times
          CasesMonth[r,q, ,d,rho_i]<-FromDaysToMoth_function((length(interval)-1),sim_data$C,nSim,month0) # number of cases per month in the q group
          peakRatio[r,q,d,rho_i]<-max(CasesMonth[r,q,c(12:22),d,rho_i])/max(CasesMonth[r,q,c(1:12),d,rho_i]) # peak ratio of the q group
          
        } #end rho_i
       } # end d
      
      save(CasesMonth,file=gsub("XX",as.character(nG),namesCases))
      save(peakRatio,file=gsub("XX",as.character(nG),namesPeaks))
      save(ArrT,file=gsub("XX",as.character(nG),namesArrT))
    } # end pop loop (index q)
    semilla<-semilla+1 # change of teh seed for the stochastic simulaitons
  } # end repite loop
  
}




