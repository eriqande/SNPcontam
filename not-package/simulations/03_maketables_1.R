# short script to produce tables from thethe MCMC simulation

library(xtable)

#### Produce a latex script for a table of the absolute difference of allele frequencies
alle_tab <- allele_table(types = out_list_01$types, data = out_list_01$allele, Lvals = out_list_01$Lvals, rhovals = out_list_01$rhovals)

### Produce a matrix for a the z-value table
z_table <- MCMC_ztable(z_df = out_list_01$z, types = out_list_01$types, rhovals = out_list_01$rhovals, Lvals = out_list_01$Lvals)

### Produce a matrix for a table showing the fraction of allele frequencies that fall into the 90% interval
a_table <- MCMC_atable(allele_df = out_list_01$allele, rhovals = out_list_01$rhovals, Lvals = out_list_01$Lvals)
