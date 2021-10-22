/* Calculating the rates derivates */
void deri(int dimension, double l[], double dl[], double t)
{/* l is the lambda in a vector way, dl is the derivate in a vector way and t the time that here we don't used. */
double vm[LIMfilas][LIMcolumnas][POBLACIONES];
int i,j,m;
double auxH, transmission;

double t0;

t0= (C.t0Month - 3.0)*30.5; /* convert monthly t0 with  t=0 as 1st of Oct, to daily t0 with t=0 as 1st of January */


VM(vm);

for (i=0; i < LIMfilas; i++)
	for (j=0; j < LIMcolumnas; j++)
	   {    
                
                transmission = C.b0*(1.0+C.delta*sin(2*3.1415*(t+t0)/365.0+C.phi));
                
                auxH = vm[i][j][0]+vm[i][j][1]+vm[i][j][2];
                       

                m=(LIMcolumnas *i+ j)*EVENTOS;
                (dl+m)[0] = C.mh * vm[i][j][0];					/* susceptibles humans mortality */
                (dl+m)[1] = transmission*vm[i][j][0]*vm[i][j][1]/auxH;          /* infection of susceptibles humans */
		(dl+m)[2] = C.mh * vm[i][j][1];                                 /* infected humans mortality*/
		(dl+m)[3] = (1.0/C.IncubationPeriod)*vm[i][j][1];               /* recovery of infected humans */
		(dl+m)[4]= C.mh * vm[i][j][3];                                  /* recovered humans mortality */
		(dl+m)[5]= C.nh*auxH;                                           /* humans birth */
                


		
		};
}

