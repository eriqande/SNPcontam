# short script to produce charts from thethe MCMC simulation

library(ggplot2)
library(grid)
library(gridExtra)

#### Produce the rho boxplots
rho_plot <- MCMC_rhoplot(rho_df = out_list_01$rho, rhovals = out_list_01$rhovals)

#### Produce the allele figure
allele_plot <- MCMC_alleleplot(allele_df = out_list_01$allele, rho = .2, loci = 20)

#### Combine and save the allele and rho figures
pdf("rho_and_allele.pdf",width=5,height=6)
grid.arrange(rho_plot,allele_plot,ncol=1)
dev.off()

#### Make the Histograms of the z-values
# the code saves the graphs into 5 separate files
# width and height denote the dementions of the pdf files
MCMC_hist(z_df = out_list_01$z, types = out_list_01$types, width = 11, height = 6)
