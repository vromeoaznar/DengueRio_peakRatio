/* here rates that varies on time */




void coeficientes(float Spark) /* Produce the coefficients that varying per day */
{
   /* double ln_sparks;
    double pop;
    pop=(double) C.Nh;
    ln_sparks=-13.2+0.00075*pop+0.85*log(pop)+(0.74-0.000143*pop)*log(Spark);
    C.spark= exp(ln_sparks)/30.5;*/
    C.spark=Spark;
    
}
