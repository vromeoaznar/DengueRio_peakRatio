/* 

Names
poblations: *) N: humans
            *) S: susceptibles humans
            *) I: infected humans
            *) R: recovered humans

rates: *) IncubationPeriod
        *) nH = mH : birth and death rates 
        *) beta = transmission rate = b0(1+delta*sin(w*t+phi)) with w=2*3.1416/365

Associated equations

dS/dt= -beta*I*S/N + nH*N - mh*S 
dI/dt= beta*I*S/N - (1/IncubationPeriod)*I - mh*I 
dR/dt= (1/IncubationPeriod)*I - mh*R
*/ 

#include "dengue.h"
#include "algoritmos.c" /* generacion de numeros aleatorios */ 

/* global variables */
common C;

double pobla[LIMfilas][LIMcolumnas][POBLACIONES];
double neventos[LIMfilas][LIMcolumnas][EVENTOS]; 
double Lbd[LIMfilas][LIMcolumnas][EVENTOS];

#include "rates_dengue.c" /* rates that are not constant in te time */
#include "f-auxiliares_dengue.c" /*coefficients and inputs rutins  */
#include "vm_dengue.c" 
#include "rk2.c"
#include "deriv_dengue.c"
//#include "reparto_dengue.c"
#include "reparto_dengue_determinista.c" // deterministic version


int main()
{
   
   
   long idum= 3480;  /*Semilla para la generacion de numeros pseudoaleatorios*/

   double t=0.0, *pointerLbd;

/* POPULATIONS
   0 susceptibles humans (S)
   1 infected humans (I)
   2 recovered humans (R) */

/* EVENTS
   0 - susceptibles humans mortality
   1 - infection of susceptibles humans
   2 - infected humans mortality
   3 - recovery of infected humans
   4 - recovered humans mortality
   5 - humans birth
  
*/

    char *poblaname[POBLACIONES], mastername[40], *eventoname[EVENTOS]; 
   char poblaN[POBLACIONES][48], eventoN[EVENTOS][48]; 
   FILE *FilePobla[POBLACIONES], *FileEventos[EVENTOS]; 
/*lectura de datos e impresion de resultados*/

   int i,l,o,Toffset; 

double *pointereventos; 
/* transitory time (in days)  */

   int Transitorio=0, Tsimulado=0, TransitorioLeido=-1, repite;
   float Dt=1./PASO; /* time step in days */

   int Plist[POBLACIONES], Elist[EVENTOS]; 
/* list of events and populations groups to print */

   //float sparkvec[5000]; /* vector with the temperature == sparkvec values*/
   //float Spark=0.0; /* variable que tiene la temperatura == spark del dia*/
   int tiempo_inicial=0; /* time in days to initialize the calculation */

   pointerLbd =Lbd[0][0]; /* pointer to go through lambdas */
   pointereventos=neventos[0][0]; /* pointer to go through events */

/* asign the pointers to the names - No elegant */
   for (i=0;i<POBLACIONES; i++) 
      poblaname[i]=poblaN[i];
   for (i=0;i<EVENTOS; i++) 
      eventoname[i]=eventoN[i];


/* reads inmput data and open files*/

   if( (i=getdata(FilePobla, FileEventos, Plist, Elist, mastername, poblaname, eventoname, &Transitorio, &Tsimulado, &idum)) < 0 ) {
      fprintf(stderr,"Error al leer %s es %d\n", DENGUEIMPUT,i);};



   TransitorioLeido=Transitorio;

/* constants independent from temperature */
#include "CONSTANTES_dengue"

/* Repetitions for statistic */
   for (repite=1; repite <= REPITE; repite++) {
      fprintf(stderr, "Repetición %d\n",repite); 
/* populations are initialized */
     pobla_iniciales(&Transitorio, &tiempo_inicial); 
     for(i=0; i< (EVENTOS) * LIMfilas * LIMcolumnas; i++)
         *(pointereventos+i)=0.0; 
	
/* transitory */
      for (l=tiempo_inicial; l< Transitorio+tiempo_inicial ; l++) {  
         //Spark=sparkvec[l-TsaveCI]; 
         //coeficientes(Spark);
	 if((l == 0) && (tiempo_inicial != 0)) save_transitorio(pobla);
	 
     for(i=0; i< (EVENTOS) * LIMfilas * LIMcolumnas; i++) 
	*(pointereventos+i)=0; 
         for(o=0; o< PASO; o++){
            t=l+(o+1.)*Dt;
	    for(i=0; i< TOTALEVENTOS; i++)
               *(pointerLbd+i)=0; /* lambdas 0 */
	    rk2(deri, Lbd[0][0], TOTALEVENTOS, t, Dt);
	   
            updatepobla_determinista();
            
	 };
      };

/* checkint that spark value is right */
    /*  if(Transitorio+tiempo_inicial-TsaveCI > 0) {
         Spark=sparkvec[Transitorio+tiempo_inicial-TsaveCI-1];  
      }
      else {
         Spark=sparkvec[0];
      }   */


      for (i=0; (i<POBLACIONES) && (Plist[i] >=0 ); i++)
         fprintf(FilePobla[i],"%d\n", repite); 
      for (i=0; (i<EVENTOS) && (Elist[i] >=0 ); i++)
  	 fprintf(FileEventos[i],"%d\n",repite);
/* file labels */
/* SIMULATION */
      if((l == 0) && (tiempo_inicial != 0) ) save_transitorio(pobla);  /*saving initial conditions*/


      Toffset= l - (TransitorioLeido+TsaveCI);
      savedata(FilePobla,FileEventos,Plist,Elist, t=(float)Toffset); /* saving output */ 

/* loop of time after saving IC */
      for (l=tiempo_inicial+Transitorio; l< -1+Tsimulado+Transitorio+tiempo_inicial ; l++) { 
	 //Spark=sparkvec[l-TsaveCI]; /* Spark values */
         //coeficientes(Spark);
         for(i=0; i< (EVENTOS) * LIMfilas * LIMcolumnas; i++) 
              *(pointereventos+i)=0.0; 
	     
	 for(o=0; o< PASO; o++) {
 	    t=l+o*Dt;
	    for(i=0; i< TOTALEVENTOS; i++)
	       *(pointerLbd+i)=0; /* lambdas 0 */
         
	    rk2(deri, Lbd[0][0], TOTALEVENTOS, t, Dt);
	    //updatepobla(&idum);
            updatepobla_determinista();
            

         }; /*end of index o*/
      
        Toffset= l - (TransitorioLeido+TsaveCI) + 1;


      savedata(FilePobla,FileEventos,Plist,Elist, t=(float)Toffset); /*we save pob, events */
            
      };  /*end time loop, index l*/      

      Transitorio= TransitorioLeido; /* Transitorio va a ser decrementado en TsaveCI al leer el archivo de CI generado */

   }/* end loop of repetitions */

   nuevasemilla(-idum); /* we save the seed. the sign minus is for ran1 of nrc and + for ran33 */

   exit(0);
}
