# If the variables below are true, then the script will not rerun any of the simulations and will use previously 
# simulated data to makes figures and tables.  If you would like to rerun the simulations change the variables
# to FALSE
sim1_data = TRUE
sim2_data = TRUE

# Must be run in the directory that contains the directories "simulations", "supplements", 
# and "manuscript"

if(!all(file.exists("simulations", "manuscript", "supplements")))  {
  stop("You must run 01_simulation_1.R in directory that includes: \"simulations\", \"manuscript\", \"supplements\"")
}

library(SNPcontam)
library(parallel)
library(ggplot2)
library(grid)
library(gridExtra)
source('~/simulations/mixed_simulation_functions.R')
source('~/simulations/simulation_functions.R')

# Single Population Simulation
#### Run Simulation for Single Population ####
if (sim1_data == FALSE){
  #### Get data set and choose Feather H Spring as an example collection of allele freqs
  tmp_data <- swfsc_chinook_baseline
  tmp_data <- tmp_data[tmp_data$Pop == "Feather_H_sp", ]
  data <- tmp_data[,-(1:4)]  # retain only the genotypes
  
  #### Run the MCMC simulations and put output in a big Rda file
  out_list_01 <- MCMC_sims(sample_data = data,
                           N = 200,
                           Lvals = c(20,60,100,200),
                           rhovals = c(0,.025,0.075,0.2,0.5), 
                           n = 100)
}

#### Make Figures and Tables Single Population MCMC ####
#### Rho and Allele Boxplots ####
  # Produce the rho boxplots
  rho_plot <- MCMC_rhoplot(rho_df = out_list_01$rho, rhovals = out_list_01$rhovals)

  # Produce the allele figure
  allele_plot <- MCMC_alleleplot(allele_df = out_list_01$allele, rho = .2, loci = 20)

  # Combine and save the allele and rho figures
  pdf("manuscript/images/rho_and_allele.pdf",width=5,height=6)
  grid.arrange(rho_plot,allele_plot,ncol=1)
  dev.off()

#### Make the Histograms of the z-values ####
  # the code saves the graphs into 5 separate files
  # width and height denote the dementions of the pdf files
  MCMC_hist(z_df = out_list_01$z, types = out_list_01$types, width = 11, height = 6, outpath = "supplements/images")

#### Make the z-table ####
  z <- out_list_01$z
  PPlim <- 0.5 # change PPlim to change the z-value limit for a contaminated sample
  z_tab <- MCMC_ztable(z = out_list_01$z, PPlim = PPlim)
  ztab_file <- "manuscript/tables/z_table_p_one_half.tex"
  cat("{\\bf (a)} Contaminated samples \n", file = ztab_file)
  cat("\\begin{center}\n", file = ztab_file, append = T)
  make_tabular(x = z_tab$TP[-1,], dd = 3, leftcolhead = "$\\rho~~~~~~~$", rightcolshead = "Number of Loci, $L$", outfile = ztab_file, append = T)
  cat("\\end{center}\n", file = ztab_file, append = T)
  cat("{\\bf (b)} Non-contaminated samples \n", file = ztab_file, append = T)
  cat("\\begin{center}\n", file = ztab_file, append = T)
  make_tabular(x = z_tab$FP, dd = 4, leftcolhead = "$\\rho~~~~~~~$", rightcolshead = "Number of Loci, $L$", outfile = ztab_file, append = T)
  cat("\\end{center}\n", file = ztab_file, append = T)

#### Make the Allele CI Table #####
  #### Produce a matrix for a table showing the fraction of allele frequencies that fall into the 90% interval
  aa <- out_list_01$allele  # give the allele df a shorter name
  # add a column to that data frame with is T/F depending on if the true value is in or out of the credible interval
  aa$inInt <- aa$alle_freq > aa$bottomint & aa$alle_freq < aa$topint
  # break data frame up on contam prob and loc number and for each, calculate fraction of CI's that overlap the true value
  CI.overlap <- aggregate(inInt ~ contam_prob + loci_number, data = aa, FUN = mean)
  # turn that into short format:
  a_table <- acast(CI.overlap, contam_prob ~ loci_number)
  # add dimname names
  names(dimnames(a_table)) <- c("rho", "Number of Loci")
  # makes the Allele CI Table
  make_tabular(a_table, leftcolhead = "$\\rho$~~~~~~~~~~", rightcolshead = "Number of Loci, $L$", outfile = "manuscript/tables/afreq_ci_coverage.tex")

# Mixture Model Simulation
#### Run Mixture Model Simulation ####
if(sim2_data == FALSE){
  # Get the data from which baseline and mixture will be chosen and the true mixing proportions
  baseline <- swfsc_chinook_baseline
  fish_pops <- ca_fishery_props # can change row indices of this matrix to use more or less sets of fishery proportions
  
  
  # Run the MCMC simulations and put output in a big Rda file
  out_list_02 <- mixed_MCMC_sims(baseline = baseline,
                                 N = 200, 
                                 p = c(0,.025, 0.075, 0.2, 0.5),
                                 fish_pops = fish_pops, 
                                 inters = 1000, 
                                 contamination = TRUE, 
                                 less.NA = 10, n = 100)
}

#### Make Figures and Tables ####
#### Rho Plot ####
mixed_MCMC_rhoplot(rho_df = out_list_02$rho_df,
                   rhovals = out_list_02$rhovals, 
                   outpath = "manuscript/images/mixed_rho.pdf", 
                   width = 5, height = 3)

#### Z Table ####
PPlim <- 0.9 # change the PPlim value to change the z value which determines a contaminated sample
mixed_MCMC_ztable(z_df = out_list_02$z_df,
                  PPlim = PPlim,
                  outpath = "manuscript/tables/mixed_z_table.tex")

#### Get Fraction of Correctly Assigned Populations and RepUnits ####
pop_info <- mixed_MCMC_population_info(u_df = out_list_02$u_df)
# returns matrix with correct_Pop and correct_Rep for each rho value and fishery

