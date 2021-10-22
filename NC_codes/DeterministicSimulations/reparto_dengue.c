/* We make the allocation of population increases */

void reparto (int i, int j, double lam[],long  paux[][LIMcolumnas][POBLACIONES], int count[], long *idum)
{int N1=0, N2=0;
double prob, q;
q=C.rho;

/*the allocation of the positive events */
			
                        /*susceptibles humans*/
			if(lam[0] > 0)
				{prob=Lbd[i][j][0]/lam[0];
				N1=binomial(prob, count[0], idum); // mueren
				N2=count[0]-N1; // se infectan
				paux[i][j][1] += N2 ; // se suma a los infectados
				neventos[i][j][0]+= N2;//N1; // sumo a los eventos de mortalidad
				neventos[i][j][1]+= binomial(q,N2,idum);//N2; // sumo a los eventos de infeccion
				}
			/* infected humans */
			if(lam[1] > 0)
				{prob=Lbd[i][j][2]/lam[1];
				N1=binomial(prob, count[1], idum); // mueren
				N2=count[1]-N1; // se recuperan
				paux[i][j][2] += N2 ; // se suman a la poblacion de recuperados
				neventos[i][j][2]+= N1; // sumo a los eventos de muertes de infectados
				neventos[i][j][3]+= N2; // sumo a los eventos de recuperacion
				}
			/* recovered humans */
			if(lam[2] > 0)
				{//prob=Lbd[i][j][8]/lam[4];
				//N1=binomial(prob, count[4], idum); // mueren
                                N1=count[2]; 
				//N2=count[4]-N1; // deberia ser cero porque no les pasa nada mas
				neventos[i][j][4]+= N1;
				}
                        /* agrego los nacimientos */
                       
                        /* de humanos */
                        N1=count[3];
                        paux[i][j][0]+=N1;
                        neventos[i][j][5]+=N1;
}

void updatepobla (long *idum)
{
  int i,j,k, count[EVENTOS]={0};
  long Etotales=0;
  long resto[LIMfilas][LIMcolumnas][POBLACIONES], restototal=0, oldresto=-1;
  long paux[LIMfilas][LIMcolumnas][POBLACIONES];
  double Ltotal, MAXl=0.0, lam[EVENTOS];

/* hago una copia de trabajo de las poblaciones */
  for (i=0; i < LIMfilas; i++)
    for (j=0; j< LIMcolumnas; j++)
      for(k=0; k<POBLACIONES ; k++)
        paux[i][j][k]=pobla[i][j][k];

  for (i=0; i< LIMfilas; i++)
    for (j=0; j< LIMcolumnas; j++) {
      Ltotal=0.0;
      for (k=0; k< (POBLACIONES); k++) {
        lam[k] = Lbd[i][j][2*k]+Lbd[i][j][2*k+1];
	resto[i][j][k]=0;
	Ltotal += lam[k];
/* lam[k] tiene la suma de los lambdas que restan a las poblacion k, tal como
esta listado en el main */
	MAXl= (MAXl >= Ltotal? MAXl : Ltotal); // que es esto?
      }
      /* le saco los nacimientos */
      lam[2] -=Lbd[i][j][5];
      Ltotal -=Lbd[i][j][5];
      Etotales= poidev(Ltotal, idum); /*total de eventos en celda, sin contar los nacimientos*/
      multinomial(0, POBLACIONES-1, Etotales, lam, count, idum);
      count[3]=poidev(Lbd[i][j][5], idum); // count[3] tiene los eventos negativos de una poblacion inexistente. Pero seran positivos para los humanos susceptibles
      
/* en count tengo los eventos negativos para cada poblacion */ // estos salen de la multinomial
      for (k=0; k< POBLACIONES; k++) {
	if (count[k] > pobla[i][j][k]) {/* exceso temporario de eventos */
	  resto[i][j][k]= count[k]-pobla[i][j][k];
	  count[k]= pobla[i][j][k];
	  paux[i][j][k]=0; /* poblacion a cero */
	  restototal += resto[i][j][k];
	}
	else /* reste no mas */ {
          paux[i][j][k] -= count[k];
	};
      };

      reparto(i, j, lam, paux, count, idum);
    }; //fin for de j e i

/* fuera del loop en i y j actualizo las poblaciones */
  for (i=0; i< LIMfilas; i++)
    for (j=0; j< LIMcolumnas; j++)
      for (k=0; k< POBLACIONES; k++)
	pobla[i][j][k]=paux[i][j][k];

  while ((restototal != oldresto) && (restototal > 0)) {
    oldresto=restototal; restototal=0;
    for (i=0; i< LIMfilas; i++)
      for (j=0; j< LIMcolumnas; j++) {
        for(k=0; k < (POBLACIONES) ; k++) {
          lam[k] = Lbd[i][j][2*k]+Lbd[i][j][2*k+1];
          count[k]=resto[i][j][k];
        };
        lam[2] -=Lbd[i][j][5]; /* nacimeintos */
     
        for(k=0; k < POBLACIONES ; k++) {
          paux[i][j][k]=pobla[i][j][k];
          if (count[k] > pobla[i][j][k]) {
            resto[i][j][k]= count[k]-pobla[i][j][k];
	    count[k]= pobla[i][j][k];
	    paux[i][j][k]=0;
	    restototal += resto[i][j][k];
          }
          else {
            paux[i][j][k] -= count[k];
	  };
        }; //fin del for k's 2

        reparto(i, j, lam, paux, count, idum);

      };// fin for j

      for (i=0; i< LIMfilas; i++)
        for (j=0; j< LIMcolumnas; j++)
	  for (k=0; k< POBLACIONES; k++)
	    pobla[i][j][k]=paux[i][j][k];
  }; /* while loop */
}//End update



