# The reversible-jump MCMC for continuous-time hidden Markov models with an unknown number of states and an unknown number of clusters 
This repository contains the R code supporting the paper by Yu Luo and and David A. Stephens (2021) ["Bayesian inference for continuous-time hidden Markov models with an unknown number of states"](https://link.springer.com/article/10.1007/s11222-021-10032-8), published in 
<em><strong>Statistics and Computing</strong></em>.

Each folder contains the code for individual simulation for: 
-	<strong>UnknownStates_Norm</strong>: RJMCMC for continuous-time hidden Markov models with an unknown number of states (Normal case, sigma=1, Example in Section 5.1 in the main paper)
-	<strong>UnknownStates_Pois</strong>: RJMCMC for continuous-time hidden Markov models with an unknown number of states (Poisson case, Example in Section 5.1 in the main paper)
-	<strong>UnknownStatesClusters_Norm</strong>: RJMCMC for continuous-time hidden Markov models with unknown numbers of states and clusters (Normal case, sigma=1, Example in Section 5.4 in the main paper)
-	<strong>UnknownStatesClusters_Pois</strong>: RJMCMC for continuous-time hidden Markov models with unknown numbers of states and clusters (Poisson case, Example in Section 5.4 in the main paper)

Please run each folder in the order: data, func, RJMCMC:
-	<strong>Data</strong>: code to generate the simulated data
-	<strong>Func</strong>: some preliminary functions, initial values and hyper parameters
-	<strong>RJMCMC</strong>: code to run the proposed RJMCMC  
