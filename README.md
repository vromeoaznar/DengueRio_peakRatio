# DengueRio_peakRatio

********** Deterministic Simulations ***************************
* The C code is a deterministic adaptacion of a code written for random simulations (poisson process). 
compile as:
gcc -lm -o dengue dengue.c or gcc -o dengue dengue.c -lm

* dengue.inp is an input file with the following structure
base_file_name_to_use_for_outputs
Number_of_populations Number_of_events
List of populations
List of events
Transitory number_of_days_of_Simulation Number_of_humans_per_unit

Variables you may change 
number_of_days_of_Simulation is a integrer number
Number_of_humans_per_unit is an float number
		
* Some parameters are given (or declared) at compile time and are defined in dengue.h (no need to change anything here to produce figures of the article )

* values for constant rates are given in CONSTANTES_dengue. Of particular interest (to produce values for figure 3 ) is C.t0Month (t0Month has units of time in months where t0Month=0 corresponding to 1st Oct)
the (Number_of_humans_per_unit,t0Month) sets used to produce Fig 3 are:
(20, 5.170492)
(45, 4.744262)
(55, 4.44918)
(119, 3.760656)
(203, 3.432787)
(400, 3.039344)
(550, 2.87541)
(675, 2.777049)
(850, 2.629508)
(1200, 2.547541)
(3000, 2.285246)
these set of values are in DeterministicValues.Rdata, as well as the corresponding peak ratio.


Obs. The C.toMonth value is provided in the file CONSTANTES_dengue. Since the deterministic runs use October as month 0 and the empirical dataset October is month 1, Fig 3B uses t0Month+1.


The code outputs are the number of events per day (X type, the types are definied in dengue.c) (output file: base_file_name_to_use-EX.dat) and number of suceptibles (base_file_name_to_use-P0.dat), infectious (base_file_name_to_use-P1.dat) and recovered (base_file_name_to_use-P2.dat) people per day per unit.


AnalysisRunsSIR_determ_C.R reads base_file_name_to_use-P1.dat and computes the monthly peak ratio (see values in DeterministicValues.Rdata)


*************** Random Simulations **********************
SIR_SparkSpatialModel_Prob.R is the code (main body) for an SIR stochastic model (demographic noise) with the different events represented as Poisson processes.  
It is implemented in R (3.6.3)with the package POMP (2.8) 

The code runs through:
 - different partitions of units according to their population (12, 25, 50 and 100 groups defined by population quantiles). The number of partitions is given by the variable nG_vector.
 - replications of the simulations. The number of replications is given by the variable Repite
 - different values of delta (this parameter is related to the amplitud of transmission seasonality; transmission_rate= b0*(1+delta*sin(w*t+phi))
 - different values of reporting rate (rho). 
obs: For possible runtime concerns, it is easy to set the values of these variables (the exploration of different values are implemented by for loops)

Runtime to simulate the 11247 of units that compound the city of Rio with 1 value for delta, 1 value for rho, and 1 replication (on regular laptop with ubuntu 16.04.7) is approximately:
4 minutes and 45 seconds for units partitioned into 12 population groups.
3 minutes and 11 seconds for units partitioned into 25 population groups.
3 minutes and 7 seconds for units partitioned into 50 population groups.
3 minutes and 8 seconds for units partitioned into 100 population groups.


Inputs: Population per unit, and sparks per unit and time (units in same population group receive the same sparks. Input sparks have units of time in months and then become daily by distributing them uniformly). These dataset are provided by the following files: 
PopMatrixNoMaxLimit_nG12.dat
PopMatrixNoMaxLimit_nG25.dat
PopMatrixNoMaxLimit_nG50.dat
PopMatrixNoMaxLimit_nG100.dat
(PopMatrixNoMaxLimit_nGXX.dat: columns correspond to population groups definied by quantiles and rows correspond to units)
nSparksMonth_nG12.dat
nSparksMonth_nG25.dat
nSparksMonth_nG50.dat
nSparksMonth_nG100.dat
(nSparksMonth_nGXX.dat: columns correspond to population groups defined by quantiles and rows correspond to month)


Outputs: Arrival time, Cases per month and peak ratio.


******* Figures *******

The Figures are produced with the file "ScriptForFigures.R"


****** the sparks ******

sparks are computed with "CompOfSparks.R". These produced datasets are used as input for the random simulations






