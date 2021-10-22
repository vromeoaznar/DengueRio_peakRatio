/* binomial.c
 * 
The binomial distribution has the form,

   prob(k) =  n!/(k!(n-k)!) *  p^k (1-p)^(n-k) for k = 0, 1, ..., n

This routine is a direct implementation of the algoritm in:

The Generation of Binomial Random Variates,  Wolfgang Hörmann,
Department of Applied Statistics and Data Processing
Wirtschaftsuniversität Wien, Preprint Series.  Preprint 1, April 1992.
http://statistik.wu-wien.ac.at/

The generation of binomial random variates
Wolfgang Hörmann
Journal of Statistical Computation and Simulation, 1563-5163, Volume 46,
Issue 1, 1993, Pages 101 – 110

Copyright 2009: H G Solari and M J Otero

 */





unsigned long binomial(double p, long n, long *idum)
{unsigned int x;
  if ((n==0) || (p==0.)) return 0;
  if (p==1.) return n;

  if (((n*p)<10.0) || (n*(1.-p))<10.) {
                  x=BINV(p,n, idum);
		  return x;
                  }
  else {
        x=BTRD(p,n, idum);
	return x;
        }

}


unsigned long BTRD(double p, long n,long *idum)
{

/*This is the algorithm from Wolfgang Hörman, Inst. f. Statistik Wirtschaftuniv. Austria */

double fc(int);         /*Stirling*/
int signo(double);      /*funcion signo con argumento no entero*/

double u=0.0; 
double v=0.0; 
int k,m,km, nm,nk;
double r,a,b,c,vr,urvr,alfa,npq,nr,t,us,h,f,ro;
int i, shift=0;;

if (p>0.5) {p=1-p;shift=1;}
m=floor((n+1)*p);
r=p/(1.-p);
nr=(n+1)*r;
npq=n*p*(1.-p);

b=1.15+2.53*sqrt(npq);
a=-0.0873+0.0248*b+0.01*p;
c=n*p+0.5;

alfa=(2.83+5.1/b)*sqrt(npq);
vr=0.92-4.2/b;
urvr=0.86*vr;


Start:
v=ran1(idum);

if(v<=urvr) 
            {u=v/vr-0.43;
            k=floor((2*a/(0.5-fabs(u))+b)*u+c);
            return shift ? n-k : k ;}

else if (v>=vr)  
        {u=ran1(idum)-0.5;}/*uniform random number in [-0.5, 0.5]*/

else {u=v/vr-0.93;
      u=signo(u)*0.5-u;
      v=ran1(idum)*vr; }/*uniform random number in [0, vr] */


/*define k from u*/
 
us=0.5-fabs(u);
k=floor((2*a/us+b)*u+c);

if ((k>=0) && (k<=n)) 
   {
   v*=alfa/(a/(us*us)+b);
   km=abs(k-m);

   if (km<=15) 
         {/*recursive evaluation of f(k)*/
         f=1;
/* DUDA */
         if (m<=k) {i=m;
                do
                {i++; f*=nr/i-r;} while(i < k);
                  }
         else if (m>k) 
                      {i=k;
		       do
                      {i++; v*=nr/i-r;} while (i< m);
                      }
         if (v<=f) return (shift ? n-k : k);
		goto Start; /* esto faltaba */
          }
   else  /*squeeze-acceptance or rejection*/
       {       v=log(v);
       ro=km/npq*(((km/3.+0.625)*km+1./6.)/npq+0.5);
       t=-km*km/(2*npq);
       if (v<(t-ro)) {return shift ? n-k : k;}
       else if (v<=(t+ro))
                          {
			  nm=n-m+1;
                          h=(m+0.5)*log((m+1)/(r*nm))+fc(m)+fc(n-m);
                          /*final test*/
                          nk=n-k+1;
			  h += (n+1)*log(((double) nm)/nk)+(k+0.5)*log((nk*r)/(k+1))-fc(k)-fc(n-k);
                          if(v<=h) return shift ? n-k : k;
                          }
	     else
			  goto Start;
       
       
       } //end else km    

   } //end if k
goto Start;

} //end BTRD




double fc(int kk)         /*Stirling*/
{

double ff=0.0,dkk=0.;
  switch (kk)
       {
       case 0: {ff=0.08106146679532726;break;}
       case 1: {ff=0.04134069595540929;break;}
       case 2: {ff=0.02767792568499834;break;}
       case 3: {ff=0.02079067210376509;break;}
       case 4: {ff=0.01664469118982119;break;}
       case 5: {ff=0.01387612882307075;break;}
       case 6: {ff=0.01189670994589177;break;}
       case 7: {ff=0.01041126526197209;break;}
       case 8: {ff=0.009255462182712733;break;}
       case 9: {ff=0.008330563433362871;break;}
       default: {dkk=kk; ff=(1./12.-(1/360.-1./1260./((dkk+1)*(dkk+1)))/((dkk+1)*(dkk+1)))/(dkk+1);break;}
       }
return ff;       
}

int signo(double kk)
{
int ff;
  if (kk<=0) {ff=-1;return ff;}
  else {ff=1;return ff;}
}


unsigned long BINV(double p, long n,long *idum)
{
/*
   If n is small and/or p is near 0 or near 1 (specifically, if
   n*min(p,1-p) < SMALL_MEAN), then a different algorithm, called
   BINV, is used which has an average runtime that scales linearly
   with n*min(p,1-p).

   This implementation by James Theiler, April 2003.

   Additional polishing for GSL coding standards by Brian Gough.  */

  int ix;                       /* return value */
  int flipped=0;
  double q,s,np;
  int BINV_CUTOFF= 110;         /* In BINV, do not permit ix too large */
  
  
  if (p > 0.5)
              {p=1.0-p;              /* work with small p */
              flipped=1;
              }

  q=1-p;
  s= p/q;
  np=n*p;


      double f0=pow(q,n);   /* f(x), starting with x=0 */

      while (1)
        {
          /* This while(1) loop will almost certainly only loop once; but
           * if u=1 to within a few epsilons of machine precision, then it
           * is possible for roundoff to prevent the main loop over ix to
           * achieve its proper value.  following the ranlib implementation,
           * we introduce a check for that situation, and when it occurs,
           * we just try again.*/

          double f=f0;
          double u=ran1(idum);

          for (ix=0; ix<=BINV_CUTOFF;++ix)
              {if (u<f) goto Finish;
              u-=f;
              /* Use recursion f(x+1) = f(x)*[(n-x)/(x+1)]*[p/(1-p)] */
              f*=s*(n-ix)/(ix+1);
              }

          /* It should be the case that the 'goto Finish' was encountered
           * before this point was ever reached.  But if we have reached
           * this point, then roundoff has prevented u from decreasing
           * all the way to zero.  This can happen only if the initial u
           * was very nearly equal to 1, which is a rare situation.  In
           * that rare situation, we just try again.*/  
        }
    

Finish:
  return (flipped) ? (n-ix) : (unsigned int)ix;
}


