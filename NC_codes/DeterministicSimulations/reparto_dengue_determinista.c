void updatepobla_determinista(){
   int i, j, k;

   double lam[EVENTOS];
   
   for(i;i<LIMfilas;i++)
      for(j;j<LIMcolumnas;j++){
           for(k;k<EVENTOS;k++){
             lam[k]=Lbd[i][j][k];
             neventos[i][j][k]+=lam[k];
           }
         pobla[i][j][0]+=lam[5]-lam[0]-lam[1];
         pobla[i][j][1]+=lam[1]-lam[2]-lam[3];
         pobla[i][j][2]+=lam[3]-lam[4];
         
       
      } 
}
