### Analysis deterministic SIR with different N and t0 (runs from C)###

# c parameters
# delta=0.99
# r0=1.15
# gamma=17
# phi=2
# t0: el t0 we use in  C

##### function to change the format of the runs of outputs of the simulation in C (dengue.c)###
ChangeFormat<-function(pobla,REPITE,Tsim,LIMfilas){
  contador<-2;
  
  if(REPITE>1){
    aux<-array(NA,dim=c(Tsim,REPITE,LIMfilas*LIMfilas),
               dimnames=list(paste0("dia ",1:Tsim),paste0("r",1:REPITE),NULL))
    for(r in 1:REPITE){
      for(dia in 1:Tsim){
        #contador<-contador+1;
        intervalo<-c((contador+1):(contador+LIMfilas))
        aux[dia,r,]<-pobla[intervalo,];
        contador<-contador+(LIMfilas+1);
      }
      contador<-contador+1;
    }
  }
  
  if(REPITE<2){
    aux<-array(NA,dim=c(Tsim,LIMfilas,LIMfilas),dimnames=list(paste0("dia ",1:Tsim),paste0("Filas",1:LIMfilas),paste0("Column",1:LIMfilas)))
    for(dia in 1:Tsim){
      #contador<-contador+1;
      intervalo<-c((contador+1):(contador+LIMfilas))
      aux[dia,,]<-pobla[intervalo,];
      contador<-contador+(LIMfilas+1);
    }
  }
  return(aux)
}


#setwd("/home/victoria/Documents/SpatialSparkDengueModel_SIR/DengueSIR_Rio_deterministic_P")
Rep<-1 # same value as variable REPITE in dengue.h
Tsim<-950 # same value as variable number_of_days_of_Simulation  in dengue.inp
LIMfilas<-2 # same value as variable LIMfilas=LIMcolumnas in dengue.h

nameFile<-"dengue_deterministic-P1.dat" 
datos<-as.matrix(read.delim(nameFile,sep= "",fill=TRUE,header=FALSE));
Inff<-ChangeFormat(datos,Rep,Tsim,LIMfilas)[,1,1]

rm(datos)
PeakRatioDeterministic<-max(Inff[200:600])/max(Inff[1:200])
print(PeakRatioDeterministic)         
   
   
  
 
  
  # #62, 93, 37.5, 27
  # nUnitsC<-LIMfilas*LIMfilas
  # DaysMonth<-c(31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31)
  # DaysMonthC<-DaysMonth[c(1:31)]
  # #names_simC<-c("sparks_popXXX_reporting0.03.dat")
  # #names_simC<-c("dengueOnlySparks_NXXX_reporting0.03-E1.dat")
  # 
  # 
  # s_input<-array(NA,dim=c(31,12))#,dimnames = list(rownames(nSparks),colnames(nSparks)))
  # ObsSimSparksCS<-array(NA,dim=c(nMonth,Rep,nUnitsC,12))
  # 
  # 
  # 
  # phi_new<-1.159567+1.6*30.5*2*3.14/365
  # tt<-c(0:700)
  # ttM<-tt+phiM
  # t0<- -18.3
  # 
  # #aux<-1+0.9*sin((ttM+t0)*2*3.14/365)
  # aux2<-1+0.9*sin((tt+t0)*2*3.14/365+phi_new)
  # plot(tt,aux2,type="l")
  # 
  # 
  # 
