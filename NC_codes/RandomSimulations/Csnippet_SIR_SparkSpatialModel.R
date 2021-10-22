

# ---- statenames ----

statenames = c(sprintf("S%d",1:L), sprintf("I%d",1:L),sprintf("R%d",1:L),"C")
acumvarnames = c("C")
obsnames = c("Y")

# ---- paramnames ----

paramnames = c(sprintf("I0%d",1:L), sprintf("R0%d",1:L), sprintf("S0%d",1:L), 
               sprintf("A%d",1:L),sprintf("B%d",1:L),"mu_H", "mu_EI", "gamma", "rho", "omega", "phi", 
               "sigma_P", "sigma_M","Lg")#, "S0", "E0","I0","R0")

# ---- rproc ----
#Process model Csnippet
rproc <- Csnippet("
                  /* variables to go through the populations */
                  double *s=&S1;
                  double *i=&I1;
                  double *r=&R1;
                  /* variables to go through the parameters A and B */
                  const double* a=&A1;
                  const double* b=&B1;
                 
                  int j, k, POBLACIONES, EVENTOS;
                  int Ll=(int)Lg; 
                  double beta;
                  int aux; 
                  
                  POBLACIONES=3;
                  EVENTOS=6; /* 2*POBLACIONES */
                  
                  double event[EVENTOS] ;//= {0};
                  double pobla[POBLACIONES];
                  
                  double rate[EVENTOS+1];
                
                  /* values for rates that do not depend of locality */
                  rate[0]=mu_H; /*birth rate, I should call it different */
                  for(k=0; k<POBLACIONES; k++)
                    rate[2*k+1]=mu_H;
                  rate[4]=gamma;
                  rate[6]=0;
                  rate[2]=0;

                  //double dW; 
                  double Nj;
                 int DailyInf;
                 int nSparks;

                  /* Here I compute things for each unit */
                  for(j=0; j<Ll; j++) {
                    DailyInf=0;
                    Nj= s[j]+i[j]+r[j];


                    pobla[0]= s[j];
                    pobla[1]= i[j];
                    pobla[2]= r[j];



                    beta=a[j]*(1+b[j]*sin(omega*(t+273)+phi)); //638 MIRAR BIEN CUANTO HAY QUE CORRER EL BETA ... Y EL SIGNO! CHEQUEADO

                    rate[2]=beta*pobla[1]/Nj;

                    
                    /* I compute the events & update pobla*/
                    event[0]=rpois(Nj*rate[0]*dt); /*the birth */

                   for(k=0; k<POBLACIONES; k++){
                      reulermultinom(2,(int) pobla[k],&rate[2*k+1],dt,&event[2*k+1]);
                      pobla[k]+=event[2*k]-event[2*k+1]-event[2*k+2];
                   }
           

                    nSparks=(int) rpois(sp_rate*dt/rho);
             
                    DailyInf= (int) event[2];
                    aux=0;
                    if(pobla[0]>=nSparks){
                      aux=nSparks;
                    }
                    if(pobla[0]<nSparks){
                      aux=pobla[0];
                    }

                    DailyInf += aux;
                    pobla[0] -= aux;
                    pobla[1] += aux;

                   /*if(pobla[0]>=0){
                     DailyInf += (int) event[2];
                     if(pobla[0]>=nSparks){
                       pobla[0] += - nSparks;
                       pobla[1] += nSparks;
                       DailyInf += nSparks;
                     }
                     if(pobla[0]<nSparks){
                       pobla[0]=0;
                       pobla[1] += pobla[0];
                       DailyInf += pobla[0];
                     }
                   }*/

                    C+=rbinom(DailyInf, rho);

               
                    /* update of the state population variables */
                    s[j]= pobla[0];
                    i[j]= pobla[1];
                    r[j]= pobla[2];

                  } /* end of j*/
                ")

# ---- rmeas ----
rmeas <- Csnippet("
                  /*double size = 1.0/sigma_M/sigma_M;
                  Y = rnbinom_mu(size,C);*/
                  Y=C;//rbinom(C,rho);
                
                  ")

rInitVic <- Csnippet("
                     double *s = &S1;
                     double *i=&I1;
                     double *r=&R1;
                     int j;
                     int Ll=(int) Lg;
                     
                     const double* s0 = &S01;
                     const double* i0 = &I01;
                     const double* r0 = &R01;

                     /* local variables */
                     for(j=0; j<Ll; j++){
                      s[j]=s0[j];
                      i[j]=i0[j];
                      r[j]=r0[j]; 
                     }
                     C=0; ")
