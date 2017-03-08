# Salmon Integrated Population Model
This site describes an integrated population model for one or more populations of Pacific salmon. The model estimates historical run sizes by year and age class, and is capable of making forecasts for 1+ years into the future. 

# Background
The model is a modification of the [ASSESSOR](https://github.com/mdscheuerell/ASSESSOR) model developed for only one population of Pacific salmon. Changes include:

1. use of Hamiltonian Monte Carlo, as implemented in the [Stan](http://mc-stan.org/) programming language, to estimate all parameters and states;
2. modeling age composition as additive log-ratios rather than as a Dirichlet distribution; and
3. treatment of harvest estimates as stochastic.

# Requirements
At a minimum, the IPM requires the following data files:

1. observed total number of adult spawners (escapement) by population and year;
2. observed age composition of adult spawners by population and year; and
3. observed total harvest by population and year.
