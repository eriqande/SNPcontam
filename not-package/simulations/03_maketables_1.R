# short script to produce tables from thethe MCMC simulation


# this must be run in the simulations directory
load("out_list_01.rda")
source("simulation_functions.R")


library(xtable)

#### Produce a latex script for a table of the absolute difference of allele frequencies
alle_tab <- allele_table(types = out_list_01$types, data = out_list_01$allele, Lvals = out_list_01$Lvals, rhovals = out_list_01$rhovals)

### Produce a matrix for a the z-value table
z_table <- MCMC_ztable(z_df = out_list_01$z, types = out_list_01$types, rhovals = out_list_01$rhovals, Lvals = out_list_01$Lvals)



#### Produce a matrix for a table showing the fraction of allele frequencies that fall into the 90% interval  ####
aa <- out_list_01$allele  # give the allele df a shorter name

# add a column to that data frame with is T/F depending on if the true value is in or out of the credible interval
aa$inInt <- aa$alle_freq > aa$bottomint & aa$alle_freq < aa$topint

# break data frame up on contam prob and loc number and for each, calculate fraction of CI's that overlap the true value
CI.overlap <- aggregate(inInt ~ contam_prob + loci_number, data = aa, FUN = mean)

# turn that into short format:
a_table <- acast(CI.overlap, contam_prob ~ loci_number)

# still to do, format for latex...