/* save the data after a day */
void savedata(FILE *FilePobla[], FILE *FileEventos[], int Plist[], int Elist[],
	float time)	
{
int err=0,i,j;

while ((Plist[err] >= 0) && (err < (POBLACIONES)) )
	{ /* type the populations */
 	fprintf(FilePobla[err], "%f \n", time);
	for (i=0; i< LIMfilas ; i++)
		{for (j=0; j< LIMcolumnas; j++)
                        fprintf(FilePobla[err], "%f\t",pobla[i][j][err]); 
                       
		fprintf(FilePobla[err], "\n"); 
		}
	err++ ; /* advances the index */
	}


err=0;  /* type the events */
while ((Elist[err] >= 0) && (err < EVENTOS) ) // TENIA B), DONDE DICE EVENTOS IBA EVENTOS+7
	{
	fprintf(FileEventos[err], "%f \n", time);
	for (i=0; i< LIMfilas ; i++)
		{for (j=0; j< LIMcolumnas; j++)
                        fprintf(FileEventos[err], "%f\t",neventos[i][j][err]);
                      
		fprintf(FileEventos[err], "\n");
		}
	err++ ; 
	} /* advances the index */
}

int getdata(FILE *FilePobla[], FILE *FileEventos[], int Plist[], int Elist[],
	char *mastername, char *poblaname[], char *eventoname[],//
	int *Transitorio, int *Tsimulado, long *idum)
{
int Npobla=0, Neventos=0;
int err=0,i;
char aux2[10]="NADA", aux1[50]="NADA";
FILE *inputfile;
inputfile=fopen(DENGUEIMPUT,"r");
if(inputfile == (FILE *) NULL)
	return -1; /* error of reading the file */
fscanf(inputfile,"%s\n%d %d",mastername,&Npobla,&Neventos);
fprintf(stderr,"Lei %s\n%d %d\n",mastername,Npobla,Neventos);

/* Open the population files */
if ((Npobla < 0) || (Npobla > POBLACIONES) )
	return err=-2;
else
	{for (i=0; i< Npobla; i++)
		fscanf(inputfile,"%d", &Plist[i]);
	} 
	for (i=Npobla; i< POBLACIONES; i++)
		Plist[i]=-1;

for (i=0; i < Npobla ; i++) 
	{
	strcpy(aux1,mastername);/* copy the master name to the auxiliary */
	
	sprintf(aux2,"-P%d.dat",Plist[i]); /* I add the indicator */
	strcat(aux1,aux2); /* I add the indicator */
	strcpy(poblaname[i],aux1); /* save the name */
	FilePobla[i]=fopen(poblaname[i],"w"); /* pass the pointer to the filepointers vector */
	};


/* the event section */
if ((Neventos < 0) || (Neventos > EVENTOS) ) // PASABA B) DONDE DICE EVENTOS IBA EVENTOS+7
	return err=-3;
else
	{for (i=0; i< Neventos; i++)
		fscanf(inputfile,"%d", (Elist+i));
	}
	for (i=Neventos; i< EVENTOS; i++) 
		Elist[i]=-1;


for (i=0; i < Neventos ; i++) 
	{
	strcpy(aux1,mastername);/* copy the master name to the auxiliary */
	
	sprintf(aux2,"-E%d.dat",Elist[i]);
	strcat(aux1, aux2); /* I add the indicator */
	strcpy(eventoname[i],aux1); /* save the name */
	FileEventos[i]=fopen(eventoname[i],"w"); /* pass the pointer to the filepointers vector*/
	};
fscanf(inputfile,"%d %d %i", Transitorio, Tsimulado, &C.Nh); 

fclose(inputfile);

if ((inputfile=fopen("semilla.inp","r") )!= NULL)
	fscanf(inputfile, "%ld", idum);
else
	{
	printf("preciso semilla ");
	scanf("%ld", idum);
	}


return err;

}


void save_transitorio(double p[][LIMcolumnas][POBLACIONES]) {
/* p is the poinr to de poblational vector */
   
   FILE *tedio;
   int i,j,k;

   if((tedio=fopen("CIniciales.dat","r")) != NULL) /* file exists */
	fprintf(stderr,"CIniciales existia\n");
	
   else {
     fprintf(stderr,"Guardo transitorio  de %d días\n", TsaveCI);
     tedio=fopen("CIniciales.dat","w");
     for (i=0;i<LIMfilas; i++) {
        for(j=0;j<LIMcolumnas;j++)
	  for(k=0; k< POBLACIONES; k++)
             fprintf(tedio,"%f ", p[i][j][k]);
    
        fprintf(tedio,"\n");
       };
  
     fclose(tedio);
   };
}

void nuevasemilla(long idum) {
   FILE *tedio;
   tedio=fopen("semilla.inp","w");
   fprintf(tedio,"%ld\n", idum);
   fclose(tedio);
}

int NumberOfVectors(int N, float Vi){
  int vectors;
  float aux;

 
  aux=Vi*N;//*pow(N*1.0,1.5);
  //aux=Vi*sqrt(N);
  vectors=(int)rint(aux);
  return vectors;
}

void pobla_iniciales(int *Transitorio, int *tiempo_inicial) {
  FILE *tedio;
  int i,j,k;

  if((tedio=fopen("CIniciales.dat","r")) == NULL) {
    fprintf(stderr,"Transitorio largo\n");
    *tiempo_inicial= TsaveCI;
    for (i=0;i<LIMfilas; i++)
      for(j=0;j < LIMcolumnas;j++) {
        pobla[i][j][0]=C.Nh*1.0;
	for(k=1; k< POBLACIONES; k++)
	  pobla[i][j][k]=0.;
	// peaks study//
        pobla[i][j][0]+= -1.;
        pobla[i][j][1] = 1.;
        /// end peaks ///	
       /* //for(k=1; k < POBLACIONES; k++)
          //pobla[i][j][k]=0;
        pobla[i][j][0]=NumberOfVectors(C.Nh, C.V0);// ACA PONER EL VALOR DE MOSQUITOS V QUE VAYAS A USAR
        pobla[i][j][1]=0;
        pobla[i][j][2]=(C.Nh-C.Inf0);
        pobla[i][j][3]=C.Inf0;
        pobla[i][j][4]=0;*/
      }; /* initialize all in zero less the susceptible mosquitoes and humans susceptibles and infected*/
    }
  else {
    *tiempo_inicial=0.;
    *Transitorio= *Transitorio+TsaveCI;
    for (i=0;i<LIMfilas; i++){
      for(j=0;j<LIMcolumnas;j++){
	for(k=0; k< POBLACIONES; k++)
          fscanf(tedio,"%lf", pobla[i][j]+k);
          
      }
    };  /* reading the initial population */
  
  }; /* reading the data of tedio */

        
  if((tedio != NULL))
    fclose(tedio); /* closing tedio */
}/* returning with the initial populations */


   
