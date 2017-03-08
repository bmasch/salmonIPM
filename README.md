# Salmon Integrated Population Model
This site describes an integrated population model for one or more populations of Pacific salmon. The model estimates historical run sizes by year and age class, and is capable of making forecasts for 1+ years into the future. 

# Background
The model is a modification of the [ASSESSOR](https://github.com/mdscheuerell/ASSESSOR) model developed for only one population of Pacific salmon. Changes include:

1. use of Hamiltonian Monte Carlo, as implemented in the [Stan](http://mc-stan.org/) programming language, to estimate all parameters and states;
2. the optional ability to model multiple populations simultaneously in a hierarchical framework;
3. allowing hatchery-origin adults to contribute to the spawning population, and distinguishing them from natural-origin adults; 
4. including removals of natural-origin adults due to hatchery broodstock collection in addition to harvest; and
5. modeling age composition as additive log-ratios rather than as a Dirichlet distribution.

# Requirements
At a minimum, the IPM requires the following data files:

1. observed total number of adult spawners (escapement) by population and year;
2. observed age composition of adult spawners by population and year (frequencies by age, based on finite samples);
3. observed proportion of total spawners that were of hatchery origin (frequencies of each rearing origin, based on finite samples); 
4. observed numbers of adults removed for broodstock by population and year; and
5. observed total harvest mortality by population and year.

Optionally, time series of environmental conditions may also be used as covariates of annual cohort productivity. In the multi-population case, these covariates are assumed to influence all populations identically. 
