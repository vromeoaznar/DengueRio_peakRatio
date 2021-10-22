############################## Fig 1 #############################

### Tot cases in the city vs time
rm(list=ls())
load("../Documents/NC_codes/Datasets/TotCases.Rdata")
par(cex=2.00)
x<-round(seq(2010,2015,1/12),1)
x<-x[-length(x)]
plot(x,TotCases,pch=1,lwd=3,xlab="Time", ylab="Cases in the city")
seasonDENV1<-c(10:21)
x_seasonDENV1<-x[seasonDENV1]
points(x_seasonDENV1,TotCases[seasonDENV1],pch=19,col=2)
seasonsDENV4<-c(22:45)
x_seasonsDENV4<-x[c(22:45)]
points(x_seasonsDENV4,TotCases[seasonsDENV4],pch=19,col=4)

par(cex=2)
legend(2013.40,max(TotCases)+500,c("DENV1", "DENV4"),pch=c(19,19),col=c(2,4), 
       x.intersp=0.9, y.intersp=0.95,adj = 0.085)#,bty = "n")
box()

##### the three incidence patterns 
rm(list=ls())

dataInterest<-read.table("../Documents/NC_codes/Datasets/CasesAggrByAR_10regions.dat")
RatioRegions<-apply(dataInterest[,c(13:24)],1,max)/apply(dataInterest[,c(1:12)],1,max)

x<-round(seq(2010,2015,1/12),1)
x<-x[-length(x)]
# Region 3.1
aux<-which((RatioRegions>0.9) & (RatioRegions<1.1))
naranjas<-colors()[grep("orange", colors())]
nNaranja<-6
par(cex=2)
plot(x[c(21:44)],dataInterest[aux,c(1:24)],pch=19,xlab="Time", ylab="Cases",col=naranjas[nNaranja])

# region 5.1 blue
aux<-which(RatioRegions==min(RatioRegions))
azules<-colors()[grep("blue", colors())]
nAzul<-24
par(cex=2)
plot(x[c(21:44)],dataInterest[aux,c(1:24)],pch=19,xlab="Time", ylab="Cases",col=azules[nAzul])

# region 2.1 green
aux<-which(RatioRegions==max(RatioRegions))
verdes<-colors()[grep("green", colors())]
nGreen<-28
par(cex=2)
plot(x[c(21:44)],dataInterest[aux,c(1:24)],pch=19,xlab="Time", ylab="Cases",col=verdes[nGreen])


############################## Fig 2 #############################

rm(list=ls())

### Peak Ratio vs Pop - aggregation by population group -
PopsMeans<-read.table("../Documents/NC_codes/Datasets/MeanPop_AllGroups.dat",header = TRUE)

peaksRatios<-array(NA,dim=dim(PopsMeans),dimnames = list(NULL,colnames(PopsMeans)))
ng<-c(12,25,50,100)
nameFile<-"../Documents/NC_codes/Datasets/CasesPerPopGroup_DENV4_nGXX.dat"

for(g in 1:length(ng)){
  CasesPerPopGroup<-read.table(gsub("XX",as.character(ng[g]),nameFile))
  peaksRatios[c(1:ng[g]),g]<-apply(CasesPerPopGroup[c(13:24),],2,max,na.rm=TRUE)/apply(CasesPerPopGroup[c(1:12),],2,max,na.rm=TRUE)
}

Mx<-max(PopsMeans,na.rm=TRUE)
mx<-min(PopsMeans,na.rm=TRUE)
My<-max(peaksRatios,na.rm=TRUE)
my<-min(peaksRatios,na.rm=TRUE)

par(cex=2)
plot(c(mx,Mx),c(my,My),type="n",xlab = "Pop density", ylab="Peak ratio")

par(cex=1.5)
points(PopsMeans[,"nG_100"],peaksRatios[,"nG_100"],lwd=3,pch=1,col=1)
points(PopsMeans[,"nG_50"],peaksRatios[,"nG_50"],lwd=3,pch=2,col=2)
par(cex=1.25)
points(PopsMeans[,"nG_25"],peaksRatios[,"nG_25"],pch=19,col=4)
points(PopsMeans[,"nG_12"],peaksRatios[,"nG_12"],pch=17,col=3)
points(PopsMeans[,"nG_12"],peaksRatios[,"nG_12"],lwd=1,pch=2,col=1)


### Peak Ratio vs Pop - aggregation by administrative region -
rm(list=ls())

nRegions<-c(10,33,160)

nameFile_Cases<-"../Documents/NC_codes/Datasets/CasesAggrByAR_XXregions.dat"
nameFile_Pop<-"../Documents/NC_codes/Datasets/PopDensity_XXAR.dat"

for(g in 1:length(nRegions)){
  dataInterest<-read.table(gsub("XX",as.character(nRegions[g]),nameFile_Cases))
  MeanPop<-read.table(gsub("XX",as.character(nRegions[g]),nameFile_Pop),header = TRUE)
  peaksRatios<-apply(dataInterest[,c(13:24)],1,max,na.rm=TRUE)/apply(dataInterest[,c(1:12)],1,max,na.rm=TRUE)
  par(cex=2)
  plot(as.numeric(MeanPop), peaksRatios,xlab="Pop density", ylab="Peak Ratio",pch=1,lwd=3)
}

################## Figure 3 - Deterministic simulations - ########################
rm(list=ls())

#DeterministicValues<-as.matrix(read.table("../Documents/NC_codes/Datasets/DeterministicValues.dat",header = TRUE))
DeterministicValues<-read.table("../Documents/NC_codes/Datasets/DeterministicValues.dat",header = TRUE)
# Arr time vs pop
MeanArrT_emp<-read.table("../Documents/NC_codes/Datasets/MeanArrTime_aggrPerPop_nG12.dat", header = TRUE)
MeanPopAllGroups<-read.table("../Documents/NC_codes/Datasets/MeanPop_AllGroups.dat",header = TRUE)
MeanPop<-MeanPopAllGroups[which(is.na(MeanPopAllGroups[,1])==FALSE),1]

par(cex=2)
plot(DeterministicValues[,"Pop"],DeterministicValues[,"C.t0Month"]+1,lwd=3,pch=2,col=4,xlab="Pop. density", ylab="Arrival time")
points(MeanPop,MeanArrT_emp[,3],pch=1,lwd=3,col=1)
lines(DeterministicValues[,"Pop"],DeterministicValues[,"C.t0Month"]+1,lwd=2,col=4)

# peak ratio vs pop
CasesPerPopGroup<-read.table("../Documents/NC_codes/Datasets/CasesPerPopGroup_DENV4_nG12.dat")
peaksRatios<-apply(CasesPerPopGroup[c(13:24),],2,max,na.rm=TRUE)/apply(CasesPerPopGroup[c(1:12),],2,max,na.rm=TRUE)

par(cex=2)
plot(DeterministicValues[,"Pop"],DeterministicValues[,"peakRatio"],lwd=3,pch=2,col=4,xlab="Pop. density", ylab="Peak Ratio")
points(MeanPop,peaksRatios,pch=1,col=1,lwd=3)
lines(DeterministicValues[,"Pop"],DeterministicValues[,"peakRatio"],col=4,lwd=2)


###################### Figure 4 ###########################################3

## Fig 4B: Box plot peak ratio from simulations and from data vs pop.
rm(list = ls())

# simulated peak ratio
load("../Documents/NC_codes/RandomSimulations/PeakRatio_incT17_r01.15_delta0.2_nQ100.Rdata")

# empirical peak ratio nG 100
CasesPerPopGroup<-read.table("../Documents/NC_codes/Datasets/CasesPerPopGroup_DENV4_nG100.dat")
RatioPeaks<-apply(CasesPerPopGroup[c(13:24),],2,max,na.rm=TRUE)/apply(CasesPerPopGroup[c(1:12),],2,max,na.rm=TRUE)
rm(CasesPerPopGroup)

# Mean Pop nG100
MeanPopAllGroups<-read.table("../Documents/NC_codes/Datasets/MeanPop_AllGroups.dat",header = TRUE)
popVector<-MeanPopAllGroups[which(is.na(MeanPopAllGroups[,4])==FALSE),4]
rm(MeanPopAllGroups)

# the plot
col_aux<-c(2,3,5)
d<-4
for(rho_i in 1:3){
  
  auxPeakRatio<-peakRatio[,,d,rho_i]
  aux1<-NULL
  aux2<-NULL
  aux3<-NULL
  for(pppp in 1:length(popVector)){
    aux1<-c(aux1,auxPeakRatio[,pppp])
    aux2<-c(aux2,c(1:nrow(auxPeakRatio)))
    aux3<-c(aux3,rep(popVector[pppp],nrow(auxPeakRatio)))
  }
  auxBox<-data.frame(y=aux1,pop=round(aux3,0),sim=aux2)

  par(cex=3)
  boxplot(y~pop, data=auxBox, ylab="Peak ratio",outline=FALSE,
          las=2,axes=F,xlab="Pop. density",lwd=2,col=col_aux[rho_i],ylim=c(0,2.6))

  box()
  axis(2,at=seq(0,3,0.5),seq(0,3,0.5),las=1)
  aux<-round(popVector,0)
  axis(1,at=c(1,seq(20,100,20)),labels=as.character(aux[c(1,seq(20,100,20))]))
  par(cex=1)
  points(seq(1,100,1),RatioPeaks,pch=19,col=1)
}

### Fig 4A: Mean arrival time - simulations and data - vs pop ###
# OBS: Arr time = 1 corresponds to October.
rm(list=ls())
# simulated arr time, nG 100
load("../Documents/NC_codes/RandomSimulations/ArrT_incT17_r01.15_delta0.2_nQ100.Rdata")
# Mean Pop nG100
MeanPopAllGroups<-read.table("../Documents/NC_codes/Datasets/MeanPop_AllGroups.dat",header = TRUE)
MeanPop<-MeanPopAllGroups[which(is.na(MeanPopAllGroups[,4])==FALSE),4]
rm(MeanPopAllGroups)
# empirical arr time, nG 100 - season 3:= 2011-2012
MeanArrTime_emp<-read.table("../Documents/NC_codes/Datasets/MeanArrTime_aggrPerPop_nG100.dat")$Season_3

# Mean arr time from the simulations for delta=0.2
nQ<-100
rho<-c(0.1, 0.3, 0.5)
d<-4

ArrT_mean_sim<-array(NA,dim=c(nQ,length(rho)),
                 dimnames = list(paste0("Pop_",as.character(round(MeanPop,0))), paste0("rho_",as.character(rho))))
ArrT_aux<-ArrT[,,1,,d,]
for(rho_i in 1:length(rho)){
  for(nq in 1:nQ){
    ArrT_mean_sim[nq,rho_i]<-mean(ArrT_aux[,,nq,rho_i],na.rm=TRUE)
  }
}
rm(ArrT, ArrT_aux)

# the plot
My<-max(MeanArrTime_emp,ArrT_mean_sim,na.rm=TRUE )
my<-min(MeanArrTime_emp,ArrT_mean_sim,na.rm=TRUE )
mx<-min(MeanPop)
Mx<-max(MeanPop)

par(cex=2)
col_aux<-c(2,3,5)
plot(c(mx,Mx),c(my,My),type="n",xlab="Pop. Density", ylab="<Arrival Time>")
points(MeanPop,MeanArrTime_emp,pch=1,lwd=2)
par(cex=0.7)
for(rho_i in 1:length(rho))
  points(MeanPop,ArrT_mean_sim[,rho_i],col=col_aux[rho_i],pch=19)
par(cex=2)
legend(2000,My,c("Data", "r 0.1", "r 0.3", "r 0.5"), pch=c(1,19,19,19),col=c(1,col_aux))

# Fig 4C: ln(empirical sparks per unit) vs ln(Total Cases)
rm(list = ls())

# Total cases
load("../Documents/NC_codes/Datasets/TotCases.Rdata")
# sparks per unit, nG 12
nSparks<-as.matrix(read.table("../Documents/NC_codes/Datasets/nSparksMonth_nG12.dat"))
# Mean Pop, nG12
MeanPopAllGroups<-read.table("../Documents/NC_codes/Datasets/MeanPop_AllGroups.dat",header = TRUE)
MeanPop<-MeanPopAllGroups[which(is.na(MeanPopAllGroups[,1])==FALSE),1]
rm(MeanPopAllGroups)

# ln of variables
ln_nSparks<-log(nSparks)
ln_nSparks[which(is.infinite(ln_nSparks)==TRUE)]<-NA
lnTotCases<-log(TotCases)

# linear fit
m<-NULL
b<-NULL

for(i in 1:ncol(ln_nSparks)){
  data_aux<-data.frame(y=ln_nSparks[,i],x=lnTotCases)
  aux<-lm(y~x,data = data_aux)
  m<-c(m,aux$coefficients[2])
  b<-c(b,aux$coefficients[1])
  rm(data_aux,aux)
}

# plot
mx<-min(lnTotCases)
Mx<-max(lnTotCases)
my<-min(ln_nSparks,na.rm=TRUE)
My<-max(ln_nSparks,na.rm=TRUE)

par(cex=2)
plot(c(mx-1.25,Mx+1.25),c(my,My),type="n",xlab="ln(Ctot)", ylab="ln(Emp Sparks/Ug)")
par(cex=1.25)
who<-c(12,10,8,6,4,2)
for(i in 1:length(who)){
  points(lnTotCases,ln_nSparks[,who[i]],col=i,pch=19)
}

for(i in 1:length(who)){
  aux_x<-c(5,10)
  aux_y<-m[who[i]]*aux_x+b[who[i]]
  lines(aux_x,aux_y,lwd=3,col=i)
}

par(cex=1.50)
legend(mx-1.45,My+0.2,c("N 37","N 119", "N 325"),col=c(6,5,4),pch=rep(19,3))
legend(8.75,-4.75,c("N 597","N 904", "N 2120"),col=c(3,2,1),pch=rep(19,3))

# Fig 4D:  plot b and m vs log pop density, nG 100
rm(list = ls())
# Total cases
load("../Documents/NC_codes/Datasets/TotCases.Rdata")
# sparks per unit, nG 100
nSparks<-as.matrix(read.table("../Documents/NC_codes/Datasets/nSparksMonth_nG100.dat"))
# Mean Pop, nG100
MeanPopAllGroups<-read.table("../Documents/NC_codes/Datasets/MeanPop_AllGroups.dat",header = TRUE)
MeanPop<-MeanPopAllGroups[which(is.na(MeanPopAllGroups[,4])==FALSE),4]
rm(MeanPopAllGroups)

# ln of variables
ln_nSparks<-log(nSparks)
ln_nSparks[which(is.infinite(ln_nSparks)==TRUE)]<-NA
lnTotCases<-log(TotCases)

# linear fit
m<-NULL
b<-NULL

for(i in 1:ncol(ln_nSparks)){
  data_aux<-data.frame(y=ln_nSparks[,i],x=lnTotCases)
  aux<-lm(y~x,data = data_aux)
  m<-c(m,aux$coefficients[2])
  b<-c(b,aux$coefficients[1])
  rm(data_aux,aux)
}

par(cex=2)
plot(log(MeanPop),m,pch=1,lwd=3,ylab="Slope (m)", xlab="ln(Pop. density)")
plot(log(MeanPop),b,pch=1,lwd=3, ylab="Intercept (b)",xlab="ln(Pop. density)")
