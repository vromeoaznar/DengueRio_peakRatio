#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define DENGUEIMPUT "dengue.inp" /* imput file  */
#define LIMfilas 2 /*number of rows in the grid */ 
#define LIMcolumnas 2  /*number of col in the grid */
#define PASO 12 /* number of steps in a day */
#define REPITE 1 /* Number of repetitions */
#define TsaveCI -0 /* Transition time (time before the simulation) */

typedef struct Common
	{
         float IncubationPeriod;
         float b0, phi, delta; 
         float mh, nh;
         float spark, rho;
	
         int Inf0, Nh;
         double t0Month;
  	} common;  /* define el type de common  for Common structure which contains the poiners to variables. These variables are used by subrutines*/
 /* some of these values are in CONSTANTES*/


#define JULIO 184         /* off set del 1 de julio a 1 de enero */ 
#define EVENTOS 6        /*total number of events per unit*/
#define POBLACIONES  3	  /* total number of populations per unit */
#define TOTALPOBLA LIMcolumnas * LIMfilas * POBLACIONES /* total number of populatins in the grid */
#define TOTALEVENTOS   LIMcolumnas * LIMfilas * EVENTOS /* total number of events in the grid */

/* funtions prototypes */
float ran1(long *);         /* uniform random number */
float gammln(float );       /* gamma distribution. we need it for the Poisson*/
int poidev(float , long *); /* Poisson distribution*/
float gasdev(long *);       /* random numbers with N(0,1) distribution*/

/* multinomial simple */
unsigned long BINV(double, long ,long *);
unsigned long BTRD(double ,long ,long *);
unsigned long binomial(double ,long ,  long *);
void multinomialS(int , int , int ,double [],int [],long *);
void multinomialR(int ,int ,int ,double [],int [],long *);



