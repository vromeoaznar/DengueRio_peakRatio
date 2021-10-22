require(maptools)
require(rgdal)
require(sp)
library(tidyverse)
load("/home/victoria/Documents/SpatialSparkDengueModel_SIR/FiguresPaper/Fig1A/DataRio.Rdata")
file_pop<-"/home/victoria/Documents/Rio_sparks/Rio_pop/Rio_poptotal_urban_grid_wgs84.shp"
DataPop<-readShapePoly(file_pop)@data

file_region<-"/home/victoria/Documents/Rio_maps/Rio_RegionsScaling/grid_bairro_cap_ra_wgs84.shp"
region<-readShapePoly(file_region)@data
aux<-duplicated(region$id)
region<-region[-which(aux==TRUE),]

regionPop<-DataPop %>% left_join(region, by='id')
regionPop<-regionPop[-which(DataPop$POP<11),]
regionPop<-regionPop[,c(5,7,8,9,10)]

regionDataRio <- DataRio %>% left_join(regionPop, by='id') #dataset with pop, cases and administrative regions

colnames(regionDataRio)<-gsub("Y","",colnames(regionDataRio))
colnames(regionDataRio)<-gsub("X","",colnames(regionDataRio))
colnames(regionDataRio)[c(32:60)]<-gsub("M","",colnames(regionDataRio)[c(32:60)])
# dataset for each AR #
# first we need to compute the max number of units in the ARs 
regionDataRio$COD_10AR<-gsub("AP ","",regionDataRio$CAP_SMS)
#namesAR<-unique(regionDataRio$COD_10AR)
namesAR<-unique(regionDataRio$BAIRRO)

nU_byAR<-NULL
for(ar in 1:length(namesAR))
  nU_byAR<-c(nU_byAR,nrow(filter(regionDataRio,BAIRRO==namesAR[ar])))
  #nU_byAR<-c(nU_byAR,nrow(filter(regionDataRio,COD_10AR==namesAR[ar])))

intDENV4<-c(22:45)
CasesByAR_DENV4<-array(NA,dim=c(max(nU_byAR),length(intDENV4),length(nU_byAR)),
                       dimnames = list(NULL,colnames(regionDataRio[,intDENV4]),namesAR))
Pop_byAR<-array(NA,dim=c(max(nU_byAR),length(nU_byAR)),
                dimnames = list(NULL,namesAR))


for(ar in 1:length(namesAR)){
  #aux<-filter(regionDataRio,COD_10AR==namesAR[ar])
  aux<-filter(regionDataRio,BAIRRO==namesAR[ar])
  CasesByAR_DENV4[c(1:nrow(aux)),,ar]<-as.matrix(aux[,intDENV4])
  Pop_byAR[c(1:nrow(aux)),ar]<-as.numeric(aux$Pop)
}
rm(aux)

# namesInt<-colnames(regionDataRio)[intDENV4]
# aggrCasesByAR<-array(NA,dim=c(length(namesAR),length(intDENV4)),
#                      dimnames = list(paste0("ARname_",namesAR),namesInt))
# for(tt in 1:length(intDENV4)){
#   aggrCasesByAR[,tt]<-apply(CasesByAR_DENV4[,tt,],2,sum,na.rm=TRUE)
# }

write.table(aggrCasesByAR,file="CasesAggrByAR_160regions.dat")
# defining pop groups -whole city- #
Pop<-regionDataRio$Pop
nQ<-12
q<-quantile(Pop,seq(0,1,1/nQ))
q[1]<-q[1]-1
HistPop<-hist(Pop,plot=FALSE,breaks = q)
mids_pop<-HistPop$mids#hist(Pop,breaks = q,plot=FALSE)$mids
breaks_pop<-q # check this
## calculating mean pop per group
MeanPop<-rep(NA,length(mids_pop))
for(gP in 1:length(mids_pop)) {
  aux<-Pop[which((Pop>breaks_pop[gP]) & (Pop<=breaks_pop[(gP+1)]))]
  MeanPop[gP]<-mean(aux)
}
# pop group per AR
nQ<-12
breaks_pop_AR<-array(NA,dim=c(length(namesAR),(nQ+1)),dimnames = list(namesAR,NULL))
mids_pop_AR<-array(NA,dim=c(length(namesAR),nQ),dimnames = list(namesAR,NULL))
for(ar in 1:length(namesAR)){
  q<-quantile(Pop_byAR[c(1:nU_byAR[ar]),ar],seq(0,1,1/nQ))
  q[1]<-q[1]-1
  HistPop_AR<-hist(Pop_byAR[c(1:nU_byAR[ar]),ar],plot=FALSE,breaks = q)
  mids_pop_AR[ar,]<-HistPop_AR$mids
  breaks_pop_AR[ar,]<-HistPop_AR$breaks
}
MeanPop_AR<-array(NA,dim=c(length(namesAR),nQ),dimnames = list(namesAR,NULL))
for(ar in 1:length(namesAR)){
  for(gP in 1:nQ){
    PopAux<-Pop_byAR[c(1:nU_byAR[ar]),ar]
    aux<-PopAux[which((PopAux>breaks_pop_AR[ar,gP]) & (PopAux<=breaks_pop_AR[ar,(gP+1)]))]
    MeanPop_AR[ar,gP]<-round(mean(aux),0)
    rm(aux,PopAux)
  }
}
rm(gP,ar)
## Now we separet by Pop group##
CasesAgg_byPop_AR<-array(NA,dim=c(length(intDENV4),nQ,length(namesAR)),
                         dimnames = list(colnames(CasesByAR_DENV4[,,1]),
                                         NULL,
                                         namesAR))

nU_byPop_byAR<-array(NA,dim=c(nQ,length(namesAR)),
                     dimnames = list(NULL,
                                     namesAR))
for(ar in 1:length(namesAR)){
  subSetPop<-Pop_byAR[c(1:nU_byAR[ar]),ar]
  for(gP in 1:nQ){
    whoPop<-which((subSetPop>breaks_pop_AR[ar,gP]) & (subSetPop<=breaks_pop_AR[ar,(gP+1)]))
    nU_byPop_byAR[gP,ar]<-length(whoPop)
    CasesAgg_byPop_AR[,gP,ar]<-apply(CasesByAR_DENV4[whoPop,,ar],2,sum,na.rm=TRUE)
    rm(whoPop)
  }
  rm(subSetPop)
}
rm(gP,ar)

peakRatio_byAR<-array(NA,dim=c(nQ,length(namesAR)),
                      dimnames = list(NULL,
                                      namesAR))
for(ar in 1:length(namesAR)){
  auxPeak2<-apply(CasesAgg_byPop_AR[c(13:24),,ar],2,max,na.rm=TRUE)
  auxPeak1<-apply(CasesAgg_byPop_AR[c(1:12),,ar],2,max,na.rm=TRUE)
  peakRatio_byAR[,ar]<-auxPeak2/auxPeak1
  rm(auxPeak1,auxPeak2)
}
peakRatio_byAR[which(is.infinite(peakRatio_byAR)==TRUE,arr.ind = TRUE)]<-NA

int_eff<-c(1,2,3,4,5,7)
int_eff<-c(2,3,7)

Mx<-max(MeanPop_AR[int_eff,])
mx<-min(MeanPop_AR[int_eff,])
my<-min(peakRatio_byAR[,int_eff],na.rm = TRUE)
My<-max(peakRatio_byAR[,int_eff],na.rm=TRUE)

Mx<-max(MeanPop_AR)
mx<-min(MeanPop_AR[int_eff,])
my<-min(peakRatio_byAR[,],na.rm = TRUE)
My<-max(peakRatio_byAR[,],na.rm=TRUE)
plot(c(mx,Mx),c(my,My),type="n",xlab="Pop",ylab="peakRatio")

for(ar in int_eff){
  points(MeanPop_AR[ar,],peakRatio_byAR[,ar],pch=19,lwd=2,col=ar)
  lines(MeanPop_AR[ar,],peakRatio_byAR[,ar],lwd=2,col=ar)
  #points(MeanPop_AR[ar,],peakRatio_byAR[,ar],pch=19,col=ar)
  #points(MeanPop,peakRatio_byAR[,ar],pch=19,col=ar)
  #lines(MeanPop,peakRatio_byAR[,ar],lwd=2,col=ar)
  #title(namesAR[ar])
}

# Now we separete by Pop group #
# with quantails with pop of city
CasesAgg_byPop_AR<-array(NA,dim=c(length(intDENV4),nQ,length(namesAR)),
                         dimnames = list(colnames(CasesByAR_DENV4[,,1]),
                                         paste0("Pop_",as.character(round(MeanPop))),
                                         namesAR))

nU_byPop_byAR<-array(NA,dim=c(nQ,length(namesAR)),
                     dimnames = list(paste0("Pop_",as.character(round(MeanPop))),
                                     namesAR))
for(ar in 1:length(namesAR)){
  subSetPop<-Pop_byAR[c(1:nU_byAR[ar]),ar]
  for(gP in 1:nQ){
    whoPop<-which((subSetPop>breaks_pop[gP]) & (subSetPop<=breaks_pop[(gP+1)]))
    nU_byPop_byAR[gP,ar]<-length(whoPop)
    CasesAgg_byPop_AR[,gP,ar]<-apply(CasesByAR_DENV4[whoPop,,ar],2,sum,na.rm=TRUE)
  }
}

peakRatio_byAR<-array(NA,dim=c(nQ,length(namesAR)),
                      dimnames = list(paste0("Pop_",as.character(round(MeanPop))),
                                      namesAR))

for(ar in 1:length(namesAR)){
  auxPeak2<-apply(CasesAgg_byPop_AR[c(13:24),,ar],2,max,na.rm=TRUE)
  auxPeak1<-apply(CasesAgg_byPop_AR[c(1:12),,ar],2,max,na.rm=TRUE)
  peakRatio_byAR[,ar]<-auxPeak2/auxPeak1
  rm(auxPeak1,auxPeak2)
}

peaksWithMeaning<-peakRatio_byAR
aux<-which(nU_byPop_byAR<100,arr.ind = TRUE)
peaksWithMeaning[aux]<-NA

Mx<-max(MeanPop)
my<-min(peaksWithMeaning,na.rm = TRUE)
My<-max(peaksWithMeaning,na.rm=TRUE)
plot(c(0,Mx),c(my,My),type="n",xlab="Pop",ylab="peakRatio")

for(ar in c(1:5,7:10)){
  points(MeanPop,peaksWithMeaning[,ar],pch=19,col=ar,lwd=2)
  #lines(MeanPop,peaksWithMeaning[,ar],lwd=2,col=ar)
  #points(MeanPop,peakRatio_byAR[,ar],pch=19,col=ar)
  #lines(MeanPop,peakRatio_byAR[,ar],lwd=2,col=ar)
  #title(namesAR[ar])
}
legend(500,My,namesAR[-6],col = c(1:10)[-6],pch=c(1:10)[-6])
