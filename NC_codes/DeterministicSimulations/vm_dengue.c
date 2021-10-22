/* Calculating the mean values ​​of the populations with border correction, include auxiliary functions. */

void	Gs(double Plus,double Minus,double *G0,double *G1,double *G2)
{/* calculo de las probabilidades de exeder por 0 1 y 2 un borde poblacional */
double z, y, yy;

if(Plus == 0){ *G0 = exp(-Minus);
	      *G1 = *G0*Minus;
	      *G2 = *G1*Minus/2.;
	      return;
	    }/* eventos que suman con tasa cero */

if(Minus == 0){*G0 = exp(-Plus);
		*G1=*G2=0.;
		return;
		} /* eventos que restan con tasa cero */

/* general case */
z=2.  * sqrt(Plus*Minus);
y= z/3.75;
yy= y*y;

if (y < 1) {
*G0=1.0+yy*(3.5156229+yy*(3.0899424+yy*(1.2067492
+yy*(0.2659732+yy*(0.360768e-1+yy*0.45813e-2)))));
*G0= *G0 * exp(-(Plus+Minus));
*G1= 2.*Minus*(0.5+yy*(0.87890594+yy*(0.51498869+yy*(0.15084934
+yy*(0.2658733e-1+yy*(0.301532e-2+yy*0.32411e-3))))));
*G1= *G1 * exp(-(Plus+Minus));
}
else {
y=1/y;
*G0=(exp(z-Plus-Minus)/sqrt(z))*(0.39894228+y*(0.1328592e-1
+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
+y*0.392377e-2))))))));
*G1=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
-y*0.420059e-2));
*G1=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
+y*(0.163801e-2+y*(-0.1031555e-1+y* *G1))));
*G1 =  sqrt(Minus/Plus)* *G1 * (exp(z-(Plus+Minus))/Plus);
}

*G2 = (- *G1 + Minus * *G0)/Plus;
return;
}

void fix(long n, double Plus, double Minus, double *Norm, double *correc,
double *G0, double *G1, double *G2)
/*calculo de la correccion por borde al valor medio devuelve Norma y correc*/
{
if((n < 0) ||  (n > 2)) {*Norm=1., *correc=0.; return;}
if ((Plus==0) && (Minus==0))
	{*correc=0.; *Norm=1; *G0=1.; *G1=*G2=0.; /* Evitar 0/0 */
	}
else
	{
Gs(Plus,Minus,G0,G1,G2);
switch (n){
     case 0: {*Norm = 1.-*G0-*G1-*G2;
	if(*Norm <= 0.)fprintf(stderr, "Case 0 Norma negativa  %g %g %g %g %g %g\n", *Norm,*G0, *G1,*G2, Plus, Minus);
	*correc = *G1+2* *G2;
        return; break;};
     case 1: {*Norm = 1- *G1- *G2;
	if(*Norm <= 0.)fprintf(stderr, "Case 1 Norma negativa %g\n", *Norm);
	*correc = *G2;
	return; break;};
     case 2: {*Norm = 1- *G2;
	if(*Norm <= 0.)fprintf(stderr, "Case 2 Norma negativa %g\n", *Norm);
	*correc = 0;
	   }
	  }
	}
}

/* This routine computes the mean values of the populations that are necessary for the rates derivate.*/

void VM(double vm[][LIMcolumnas][POBLACIONES])
/* Lambdas, pobla y nevents son globales. Lbd los Lambdas es global */
{
int i,j,k,im,ip,jm,jp;
long n;
double LbdP, LbdM, G0, G1, G2, Norm, correc, correc2, hold;

for (i=0 ; i< LIMfilas; i++)
	for (j=0; j < LIMcolumnas; j++)
		{

/* fix retorna la correccion por borde de poblacion */



/* Suceptibles humans */
        n=(int)pobla[i][j][0]; 
        // TENIA 
	//n=pobla[i][j][0];
	LbdP=Lbd[i][j][5]; // human birth
	LbdM=Lbd[i][j][0]+Lbd[i][j][1]; // mortality S and infection (both in and out)
	fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
	//vm[i][j][0]= n+(LbdP-LbdM + correc)/Norm;
         vm[i][j][0]= pobla[i][j][0]+LbdP-LbdM;

/* Infected humans */
        n=(int)pobla[i][j][1];
        //TENIA
	//n=pobla[i][j][1];
	LbdP=Lbd[i][j][1]; // infection (both in and out)
	LbdM=Lbd[i][j][2]+Lbd[i][j][3]; // mortality I and recovery
	fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
	hold = LbdP-LbdM;
	//vm[i][j][1]= n+(LbdP-LbdM + correc)/Norm;
        vm[i][j][1]= pobla[i][j][1]+LbdP-LbdM;

/* Recovered humans */
        n=(int)pobla[i][j][2];
        //TENIA
	//n=pobla[i][j][2];
	LbdP=Lbd[i][j][3]; // recovery
	LbdM=Lbd[i][j][4]; // mortality R
	fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
	vm[i][j][2]= pobla[i][j][2]+LbdP-LbdM;

}; // fin del for i j
} // fin de la funcion
