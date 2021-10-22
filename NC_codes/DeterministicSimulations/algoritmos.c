float randC(long *dseed)
{
/*	*dseed = random(); */
	*dseed = rand();
         return ((float) *dseed) / RAND_MAX ;
}


#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software V)&t. */

float gammln(float xx)

{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;

  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

/* Distribucion Poisson */

int poidev(float xm, long *idum)
{

   float gammln(float xm);
   static float sq,alxm,g,oldm=(-1.0);
   float em,t,y;



   if (xm < 12.0) {
      if (xm != oldm) {
         oldm=xm;
         g=exp(-xm);
      }
      em = -1;
      t=1.0;
      do {

         ++em;
         t *= ran1(idum);
         } while (t > g);
   } else {
      if (xm != oldm) {
          oldm=xm;
          sq=sqrt(2.0*xm);
          alxm=log(xm);
          g=xm*alxm-gammln(xm+1.0);

      }
      do {
         do {
            y=tan(M_PI*ran1(idum));

            em=sq*y+xm;
         } while (em < 0.0);
         em=floor(em);



         t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
      } while (ran1(idum) > t);
   }
   return em;
}
/* generacion de numeros con distribucion normal N(0,1) */

float gasdev(long *idum)

{
float ran1(long *idum);
static int iset=0;
static float gset;
float fac,rsq,v1,v2;

if (*idum < 0) iset=0;
if (iset == 0) {

   do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
      } while (rsq >= 1.0 || rsq == 0.0);
      fac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;

      iset=1;
      return v2*fac;
      } else {
        iset=0;
        return gset;


   }
}



/*Binomial del numerical recipes*/

float bnldev(float pp, int n, long *idum)
                                                                       /*Returns as a floating-point number an integer value that is a random deviate drawn from
                                                                          a binomial distribution of n trials each of probability pp, using ran1(idum) as a source of
                                                                          uniform random deviates.*/
{
     float gammln(float xx);
     float ran1(long *idum);
     int j;
     static int nold=(-1);
     float am,em,g,angle,p,bnl,sq,t,y;
     static float pold=(-1.0),pc,plog,pclog,en,oldg;
     p=(pp <= 0.5 ? pp : 1.0-pp);
                                                                                          /*The binomial distribution is invariant under changing pp to 1-pp, if we also change the
                                                                                           answer to n minus itself; we'll remember to do this below.
                                                                                           This is the mean of the deviate to be produced.*/
     am=n*p;
                                                                                            /*Use the direct method while n is not too large.*/
     if (n < 25) {
                                                                                               /*This can require up to 25 calls to ran1.*/
         bnl=0.0;
         for (j=1;j<=n;j++)
             if (ran1(idum) < p) ++bnl;
                                                                                             /*If fewer than one event is expected out of 25*/
     } else if (am < 1.0) {
                                                                                                /*or more trials, then the distribution is quite*/
         g=exp(-am);
                                                                                                   /* accurately Poisson. Use direct Poisson method.*/
         t=1.0;
         for (j=0;j<=n;j++) {
             t *= ran1(idum);
             if (t < g) break;
         }
         bnl=(j <= n ? j : n);
                                                                                                 /*Use the rejection method.*/
     } else {

                                                                                                  /*If n has changed, then compute useful quantities.*/
      if (n != nold) {
                                          
          en=n;
          oldg=gammln(en+1.0);
          nold=n;
                                                                                                   /*If p has changed, then compute useful quantities*/
      } if (p != pold) {
                                          
          pc=1.0-p;
          plog=log(p);
          pclog=log(pc);
          pold=p;
      }
                                                                                                      /*The following code should by now seem familiar:*/
      sq=sqrt(2.0*am*pc);
                                                                                                       /*rejection method with a Lorentzian comparison function*/
      do {
                                          
          do {
              angle=M_PI*ran1(idum);
              y=tan(angle);
              em=sq*y+am;
                                                                                                              /*Reject.*/
          } while (em < 0.0 || em >= (en+1.0));
                                                                                                               /*Trick for integer-valued distribution.*/
          em=floor(em);
          t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
              -gammln(en-em+1.0)+em*plog+(en-em)*pclog);
                                                                                                               /* Reject. This happens about 1.5 times per deviate, on average*/
      } while (ran1(idum) > t);
                                         
      bnl=em;
  }
                                                                                                                      /*Remember to undo the symmetry transformation*/
  if (p != pp) bnl=n-bnl;
                                          
  return bnl;
}


/* binomial.c
 * 
The binomial distribution has the form,

   prob(k) =  n!/(k!(n-k)!) *  p^k (1-p)^(n-k) for k = 0, 1, ..., n


 */


float ran1(long *);



#include "binomial.c"

void multinomialS(int is0, int ie0, int Ntiradas,double lam[],int count[],long *idum)
 {                                                                                                    /* Ntiradas numero de dados*/
                                                                                                        /*lam[ie0-is0+1]  eventos*/
                                                                                                       /*idum es la semilla del ran1 de los numeros aleatorios*/
int i, is, ie, k, event;
double au1,p[30],sumalam, padd[30];

if (Ntiradas <= 0) return;

event=ie0-is0+1;

if (event == 1) 
   {
     count[is0]=Ntiradas;             /* o bien count[is0]=count[is0]+Ntiradas*/
     return;
   }

sumalam=0;

for (i=0; i < event; i++)
      {
sumalam=lam[i+is0]+sumalam;	                                     /*sumo los valores de los eventos para luego normalizar*/
      }

for (i=0; i < event; i++)
{
p[i]=lam[i+is0]/sumalam;                                                  /*p=probabilidades, eventos normailzados por la sumalam*/     
}

au1=0;
padd[0]=p[0];
	
for (i=1; i < event; i++)
      {
         padd[i] = p[i] + padd[i-1];                                       /*tasas normalizadas acumuladas*/
      }
	                                     
	
for (k=0; k < Ntiradas; k++) 
    {
        au1=ran1(idum);
	is=0;
	ie=event;
	i=ie/2;

	do
	         {
                    if (au1 <= padd[i-1])
	               {
                          ie=i;
		  	  i=(ie+is)/2;
		        }
	 	    else
		       {
                          is=i;
		  	   i=(is+ie+1)/2;
		       }
                 }
	while ((i < ie) && (i > is));

        count[is+is0] = count[is+is0]+1;
		
     }
return;

 }


void multinomialR(int is0,int ie0,int Ntiradas,double lam[],int count[],long *idum)
 {                                                                                                                                     
	
int i, k;
double prob, sumalam, p[30];                                                                                  /*tasas normalizadas acumuladas*/ 
int N1;                                                                                        /*numero de exitos de la binomial*/
int N2;                                                                                        /*numero de fracasos de la binomial*/

prob=0;                                                       
	
for (k=is0; k<= ie0; k++)
      {
        count[k]=0;                                                   /*inicializo cada contador en cero*/  
      }

if (Ntiradas <= 0) return;

/*Si Ntiradas es menor que XX paso a la multinomialS	*/

if (Ntiradas <= 30 && Ntiradas >0)                                              
  {
     multinomialS(is0,ie0,Ntiradas,lam,count,idum);
     return;
  }

      
if (is0==ie0)                                 
   {
     count[is0]=Ntiradas;
     return;
    }

sumalam=0;

for (k=is0; k <= ie0; k++)
      {
sumalam=lam[k]+sumalam;	                                     /*sumo los valores de los eventos para luego normalizar*/
      }

for (k=is0; k <= ie0; k++)
{
p[k]=lam[k]/sumalam;                                                  /*p=probabilidades, eventos normalizados por la sumalam*/     
}

   
i=(is0+ie0)/2;                                                         /*Determino la posicion media del intervalo*/

        
for (k=is0; k <= i; k++)
      {
         prob = p[k] + prob;                                   /*Sumo las tasas normalizadas desde is hasta i, prob*/
      }
	N1= binomial(prob,Ntiradas,idum);           /*calculo el N1 correpondiente a la primer mitad del intervalo y N2 por diferencia*/
	N2=Ntiradas-N1;
        multinomialR(is0,i,N1,lam,count,idum);
        multinomialR(i+1,ie0,N2,lam,count,idum);
return;
 }


/* lam tiene los pesos no normalizados */

void multinomial(int is0,int ie0,int Ntiradas,double lam[],int count[],long *idum)
 {                                                                                                                                     
	
int i, k;
double prob, sumalam, p[30];                                                                                  /*tasas normalizadas acumuladas*/ 
int N1;                                                                                        /*numero de exitos de la binomial*/
int N2;                                                                                        /*numero de fracasos de la binomial*/

prob=0;                                                       
	
for (k=is0; k<= ie0; k++)
      {
        count[k]=0;                                                   /*inicializo cada contador en cero*/  
      }

if (Ntiradas <= 0) return;

/*Si Ntiradas es menor que XX paso a la multinomialS	*/

if (Ntiradas <= 20 && Ntiradas >0)                                              
  {
     multinomialS(is0,ie0,Ntiradas,lam,count,idum);
     return;
  }

      
if (is0==ie0)                                 
   {
     count[is0]=Ntiradas;
     return;
    }

sumalam=0;

for (k=is0; k <= ie0; k++)
      {
sumalam=lam[k]+sumalam;	                                     /*sumo los valores de los eventos para luego normalizar*/
      }

for (k=is0; k <= ie0; k++)
{
p[k]=lam[k]/sumalam;                                                  /*p=probabilidades, eventos normalizados por la sumalam*/     
}

   
i=(is0+ie0)/2;                                                         /*Determino la posicion media del intervalo*/

        
for (k=is0; k <= i; k++)
      {
         prob = p[k] + prob;                                   /*Sumo las tasas normalizadas desde is hasta i, prob*/
      }

	N1= binomial(prob,Ntiradas,idum);           /*calculo el N1 correpondiente a la primer mitad del intervalo y N2 por diferencia*/
	N2=Ntiradas-N1;
        multinomial(is0,i,N1,lam,count,idum);
        multinomial(i+1,ie0,N2,lam,count,idum);
return;

 }

