

# short script to run the simulations assessing the mixed MCMC method for
# detecting contaminated samples.

# Must be run in the directory that contains the directories "simulations", "supplements", 
# and "manuscript"

if(!all(file.exists("simulations", "manuscript", "supplements")))  {
  stop("You must run 01_simulation_1.R in directory that includes: \"simulations\", \"manuscript\", \"supplements\"")
}

library(SNPcontam)
library(parallel)
library(gtools)

source("simulations/mixed_simulation_functions.R")


#### Get the data from which baseline and mixture will be chosen and the true mixing proportions ####
baseline <- swfsc_chinook_baseline
load("data/ca_fishery_props.rda")
fish_pops <- ca_fishery_props


#### Run the MCMC simulations and put output in a big Rda file ####
set.seed(10)  # set seed for reproducibility
out_list_02 <- mixed_MCMC_sims(
                                baseline = baseline,
                                N = 200, 
                                p = c(0,.025, 0.075, 0.2, 0.5), 
                                fish_pops = fish_pops, 
                                inters = 1000, 
                                contamination = TRUE, 
                                less.NA = 10, 
                                n = 20,
                                MAX_CORES = 20)

save(out_list_02, file = "simulations/out_list_02.rda", compress = "xz")

