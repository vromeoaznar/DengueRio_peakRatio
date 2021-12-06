# DengueRio_peakRatio

********** Deterministic Simulations ***************************
* The C code is a deterministic adaptation of an implementation for random simulations (based on Poisson processes). 

compile the code as:
gcc -lm -o dengue dengue.c or gcc -o dengue dengue.c -lm

* dengue.inp is an input file with the following structure
base_file_name_to_use_for_outputs
Number_of_populations Number_of_events
List of populations
List of events
Transitory number_of_days_of_Simulation Number_of_humans_per_unit

Variables you may change in this file
number_of_days_of_Simulation is an integer number
Number_of_humans_per_unit is a float number
		
* Some parameters are given (or declared) at compile time and are defined in dengue.h (no need to change anything here to produce the figures of the article )

* values for constant rates are given in CONSTANTES_dengue. Of particular interest (to produce values for figure 3 ) is C.t0Month (C.t0Month has units of time in months where C.t0Month=0 corresponds to Oct 1st)
The sets of (Number_of_humans_per_unit, t0Month) used to produce Fig 3 are:

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

these sets and their corresponding peak ratio (i.e. the value of the peak ratio obtained for each of the sets of (Number_of_humans_per_unit, t0Month) shown above) are found in the DeterministicValues.Rdata file. 


Observation. The C.toMonth value is provided in the file CONSTANTES_dengue. Since the deterministic runs use October as month 0, and in the empirical dataset October is month 1, Fig 3B uses t0Month+1.


The output files provide the number of events per day (for X type, with the types defined in dengue.c) base_file_name_to_use-EX.dat, as well as the number of susceptible (base_file_name_to_use-P0.dat), infectious (base_file_name_to_use-P1.dat) and recovered (base_file_name_to_use-P2.dat) individuals per day per unit.


AnalysisRunsSIR_determ_C.R reads base_file_name_to_use-P1.dat and computes the monthly peak ratio (see values in DeterministicValues.Rdata)

*************** Random Simulations **********************

SIR_SparkSpatialModel_Prob.R is the code (main body) for an SIR stochastic model (with demographic noise) in which the different events are represented as Poisson processes.  
It is implemented in R (3.6.3) with the package Pomp (2.8). That is, the model is written as a Pomp object, and simulated with this package.

The code runs through:
 - different partitions of units according to their population (12, 25, 50 and 100 groups defined by population quantiles), where the number of partitions is given by the variable nG_vector.
 - replications of the simulations, with the number of replications is given by the variable Repite
 - different values of delta (this parameter is related to the amplitude of transmission seasonality, where the transmission_rate= b0*(1+delta*sin(w*t+phi))
 - different values of the reporting rate (rho). 

obs: For possible concerns regarding runtime, it is easy to set the values of these variables and check it. (The simulations for the different values are implemented via for loops)

Runtime to simulate the 11247 of units that constitute the city of Rio de Janeiro with 1 value for delta, 1 value for rho, and 1 replication (on a regular laptop with ubuntu 16.04.7) is approximately:
4 minutes and 45 seconds for units partitioned into 12 population groups.
3 minutes and 11 seconds for units partitioned into 25 population groups.
3 minutes and 7 seconds for units partitioned into 50 population groups.
3 minutes and 8 seconds for units partitioned into 100 population groups.


Inputs: Population per unit, and sparks per unit and time (units in the same population group receive the same number of sparks. Input sparks rate have units of time in months, but are implemented daily by distributing them uniformly over time). These datasets are provided by the following files: 
PopMatrixNoMaxLimit_nG12.dat
PopMatrixNoMaxLimit_nG25.dat
PopMatrixNoMaxLimit_nG50.dat
PopMatrixNoMaxLimit_nG100.dat
(PopMatrixNoMaxLimit_nGXX.dat: the columns correspond to population groups defined by quantiles and the rows correspond to units)
nSparksMonth_nG12.dat
nSparksMonth_nG25.dat
nSparksMonth_nG50.dat
nSparksMonth_nG100.dat
(nSparksMonth_nGXX.dat: the columns correspond to population groups defined by quantiles and the rows correspond to month)


Outputs: Arrival time, cases per month and peak ratio.


******* Figures *******

The Figures are produced with the file "ScriptForFigures.R"


****** the sparks ******

Sparks are computed with "CompOfSparks.R" from the observed presence and absence of case in the unit, following the approach described in the paper. The resulting datasets are used as input for the random simulations.
