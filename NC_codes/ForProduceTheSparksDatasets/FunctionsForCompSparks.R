##################### auxiliar functions for compute the sparks #######################

localOutbreaks_InitAndFinalTimes<-function(nUnits,BinaryCases,nMonth){
  nMonthCases<-apply(BinaryCases,1,sum)
  tiempos<-array(NA,dim=c(nUnits,16,2),dimnames = list(NULL,NULL,c("to","tf")))
  for(n in 1:nUnits){
    
    if(nMonthCases[n]>0){
      t_termina<-NULL
      to_guardo<-NULL
      tt<-1
      
      while(tt<(nMonth-1)){
        aux<-which(BinaryCases[n,c(tt:nMonth)]>0)
        if(length(aux)>0){
          ti<-min(aux)+tt-1
          to_guardo<-c(to_guardo,ti)
          tt<-ti
          new_season<-FALSE
          
          while(new_season==FALSE){
            tt<-tt+1
            if(tt<(nMonth)){
              if((BinaryCases[n,tt]+BinaryCases[n,(tt+1)])==0){
                ##if((BinaryCases[n,tt]+BinaryCases[n,(tt+1)])<2){
                t_termina<-c(t_termina,tt)
                new_season<-TRUE
              }
            }
            else
              new_season<-TRUE
          }
        }
        else
          tt<-nMonth
      }
      #correccion#
      if(to_guardo[length(to_guardo)]==nMonth)
        t_termina<-c(t_termina,(nMonth+1))
      if(length(to_guardo)>length(t_termina))
        t_termina<-c(t_termina,nMonth)
      
      tiempos[n,c(1:length(to_guardo)),"to"]<-to_guardo
      if(is.null(t_termina)==TRUE)
        t_termina<-nMonth
      tiempos[n,c(1:length(t_termina)),"tf"]<-t_termina
    }
  }
  
  return(tiempos)
}

########################### computing n positive untis per time ##########################
## beacause I am assuming that a local time serie with : cases, cases, 0, cases == cases,cases,cases,cases
## I will use the follow matrix to compute the number of positve units ##

missedPositiveMonth_function<-function(tiempos,MatrixCases){
  
  
  binaryCases<-array(0,dim=dim(MatrixCases),dimnames =list(NULL,colnames(MatrixCases)) )
  for(u in 1:nrow(MatrixCases)){
    aux<-which(is.na(tiempos[u,,1])==FALSE)
    if(length(aux)>0){
      interval<-tiempos[u,aux,]
      if(length(aux)>1){
        interval[,2]<-interval[,2]-1
        for(i in 1:length(aux))
          binaryCases[u,c(interval[i,1]:interval[i,2])]<-rep(1,(interval[i,2]-interval[i,1]+1))
      }
      if(length(aux)==1){
        interval[2]<-interval[2]-1
        binaryCases[u,c(interval[1]:interval[2])]<-rep(1,(interval[2]-interval[1]+1))
      }
    }
  }
  return(binaryCases)
  
}
#########      to compute the number of extiontios from the data      #########
nDynamics_function<-function(nTotMonthes, Pop, breaksPop, midsPop, tiempos, MatrixCases){
  nPopGroups<-length(breaksPop)-1
  nExtintionsByPop_output<-array(NA,dim=c(nTotMonthes,nPopGroups),
                                 dimnames=list(paste0("month=",c(1:nTotMonthes)),paste0("pop=",midsPop)))
  nPositiveUnits_output<-nExtintionsByPop_output
  nSparks_output<-nExtintionsByPop_output
  nUnitsPerPop_output<-rep(NA,length(mids_pop))
  TotPop_perGroup_output<-rep(NA,length(mids_pop))
  for(nI in 1:nPopGroups){
    aux<-which((Pop>breaksPop[nI]) & (breaksPop[(nI+1)]>=Pop))
    nUnitsPerPop_output[nI]<-length(aux)
    TotPop_perGroup_output[nI]<-sum(Pop[aux])
    
    aux1<-tiempos[aux,,"tf"]
    for(tt in 1:nTotMonthes){
      aux2<-MatrixCases[,tt]#binaryCases[,tt]#MatrixCases[,tt]
      nPositiveUnits_output[tt,nI]<-length(which(aux2[aux]>0))
      #if(tt<nTotMonthes)
      nExtintionsByPop_output[tt,nI]<-nrow(which(aux1==tt,arr.ind = TRUE))
    }
    for(tt in 2:(nTotMonthes-1)){
      nSparks_output[tt,nI]<-(nPositiveUnits_output[tt,nI]-nPositiveUnits_output[(tt-1),nI]+
                                nExtintionsByPop_output[tt,nI])/(length(aux)-nPositiveUnits_output[(tt-1),nI])
    }
  }
  
  list(nSparks=nSparks_output,nExtintionsByPop=nExtintionsByPop_output, nPositiveUnits=nPositiveUnits_output,
       nUnitsPerPop=nUnitsPerPop_output, TotPop_perGroup=TotPop_perGroup_output)
}
