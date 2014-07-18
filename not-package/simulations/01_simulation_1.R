

# short script to run the simulations assessing the MCMC method for
# detecting contaminated samples.

# Must be run in the directory that contains the directories "simulations", "supplements", 
# and "manuscript"

if(!all(file.exists("simulations", "manuscript", "supplements")))  {
  stop("You must run 01_simulation_1.R in directory that includes: \"simulations\", \"manuscript\", \"supplements\"")
}

library(SNPcontam)
library(fullsniplings)
library(parallel)

source("simulations/simulation_functions.R")

#### Get data set and choose Feather H Spring as an example collection of allele freqs ####
tmp_data <- swfsc_chinook_baseline
tmp_data <- tmp_data[tmp_data$Pop == "Feather_H_sp", ]
data <- tmp_data[,-(1:4)]  # retain only the genotypes



#### Run the MCMC simulations and put output in a big Rda file
set.seed(10)  # Set Seed For reproducibility
out_list_01 <- MCMC_sims(
                        sample_data = data,
                        N = 200,
                        Lvals = c(20,60,100,200),
                        rhovals = c(0,.025,0.075,0.2,0.5), 
                        n = 100 
                        )
save(out_list_01, file = "simulations/out_list_01.rda", compress = "xz")

