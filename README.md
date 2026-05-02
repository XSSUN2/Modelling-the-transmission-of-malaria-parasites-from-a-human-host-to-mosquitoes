The file Models.R contains the code for the mosquito feeding model and the within-human parasite dynamics model. The mosquito feeding model defined in Models.R is required for model fitting, and both models are required for simulation of asymptomatic infections and PRCC analysis, so Models.R should be run or sourced before running other codes.

Feeding assay data.csv provides the direct feeding assay data used in model fitting, including the male and female gametocytemia (measured by PfMGET and Pfs25 qRT-PCR respectively) and the proportion of mosquitoes developing oocysts.

The file ABC-SMC_vector.R provides the code for model fitting, parameter estimation and the visualization.

The file Simulation of asymptomatic infections.R provides the code for characterizing the probability of human-to-mosquito transmission of asymptomatic infections.

The file Sensitivityanalysis.R provides the code for the PRCC analysis, which identifying the the dominant host factors determining the probability of human-to-mosquito transmission.

The files Uninfected_mosquito.csv and Infected_mosquito.csv provides the mosquito survival data.
The file Death rate.R provides estimates of mosquito death rates by fitting a Gompertz function to mosquito survival data.






