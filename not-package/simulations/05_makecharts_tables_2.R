# short script to produce charts from the mixed MCMC simulation

# this script should be run from the directory that includes the directories "simulations", "supplements", 
# and "manuscript"
if(!all(file.exists("simulations", "manuscript", "supplements")))  {
  stop("You must run 01_simulation_1.R in directory that includes: \"simulations\", \"manuscript\", \"supplements\"")
}

load("simulations/out_list_02.rda")
source("simulations/mixed_simulation_functions.R")

library(ggplot2)

#### Rho Plot ####
message("     making manuscript/images/mixed_rho.pdf")
mixed_MCMC_rhoplot(rho_df = out_list_02$rho_df,
                   rhovals = out_list_02$rhovals, 
                   outpath = "manuscript/images/mixed_rho.pdf", 
                   width = 5, height = 3)

#### Z Table ####
message("     making manuscript/tables/mixed_z_table.tex")
mixed_MCMC_ztable(z_df = out_list_02$z_df,
                  PPlim = 0.9,
                  outpath = "manuscript/tables/mixed_z_table.tex")

#### Get Fraction of Correctly Assigned Populations and RepUnits ####
message("Skipping the fraction of correctly assigned fish to population and RepUnit")
message(" because that took too long...(eric commented it out)")
## Eric commented this out because it seems to take too long / hang on the large simulation output
# pop_info <- mixed_MCMC_population_info(u_df = out_list_02$u_df)