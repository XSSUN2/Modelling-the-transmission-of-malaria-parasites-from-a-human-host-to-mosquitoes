The file `Models.R` contains the code for the mosquito feeding model and the within-human parasite dynamics model. The mosquito feeding model defined in `Models.R` is required for model fitting, and both models are required for simulations of asymptomatic infections and PRCC analysis. Therefore, `Models.R` should be sourced before running other scripts (provided in the code files).

The file `Feeding assay data.csv` contains the direct feeding assay data used for model fitting, including male and female gametocytemia (measured by PfMGET and Pfs25 qRT-PCR, respectively) and the proportion of mosquitoes developing oocysts.

The file `ABC-SMC_vector.R` provides the code for model fitting, parameter estimation, and visualization using ABC-SMC.

The file `Simulation of asymptomatic infections.R` provides the code for characterizing the probability of human-to-mosquito transmission from asymptomatic infections.

The file `Sensitivityanalysis.R` provides the code for the PRCC analysis, identifying the dominant host factors determining the probability of human-to-mosquito transmission.

The estimation of mosquito death rates is provided in the `supplementary` folder.

This code was developed in R (version ≥ 4.2.0). All analyses were performed on a MacBook Pro with an Apple M1 Pro chip (8-core CPU) and 16 GB RAM. The model fitting and PRCC analysis each typically take approximately one day using parallel computation across 8 cores.



