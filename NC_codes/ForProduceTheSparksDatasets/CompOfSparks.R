# the sparks per unit per month computed here are the input of the
# random simulations ... file: "SIR_SparkSpatialModel_Prob.R"

rm(list=ls())
require("FunctionsForCompSparks.R")

datasetUnit<-read.csv("Dataset_id_Pop_NegPos.csv",header = TRUE)
MatrixCases<-datasetUnit[,c(3:ncol(datasetUnit))] # positive and negative units in time 
Pop<-datasetUnit[,"Pop"]

nUnits<-nrow(datasetUnit)
nMonth<-ncol(MatrixCases)#31 #29
tiempos<-localOutbreaks_InitAndFinalTimes(nUnits,MatrixCases,nMonth)


### Pop binned in quantiles ###
nQ<-12 # number of population groups # 12, 25, 50 and 100
q<-quantile(Pop,seq(0,1,1/nQ))
q[1]<-q[1]-1
HistPop<-hist(Pop,plot=FALSE,breaks = q)
mids_pop_quantile<-HistPop$mids
breaks_pop<-q # check this
mids_pop<-mids_pop_quantile


nDynamics<-nDynamics_function(nMonth, Pop, breaks_pop, mids_pop, tiempos, 
                              missedPositiveMonth_function(tiempos,MatrixCases))

write.table(nDynamics$nSparks,file=gsub("XX",as.character(nQ),"nSparksMonth_nGXX.dat"))

