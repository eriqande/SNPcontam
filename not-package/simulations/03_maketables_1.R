# short script to produce tables from thethe MCMC simulation


# this script should be run from the directory that includes the directories "simulations", "supplements", 
# and "manuscript"
if(!all(file.exists("simulations", "manuscript", "supplements")))  {
  stop("You must run 01_simulation_1.R in directory that includes: \"simulations\", \"manuscript\", \"supplements\"")
}

load("simulations/out_list_01.rda")
source("simulations/simulation_functions.R")


library(xtable)
library(reshape2)

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

# add dimname names
names(dimnames(a_table)) <- c("rho", "Number of Loci")

# here is a function to make output for tabular so we can customize more.  Incomplete.
# Almost have it figured out, but....later...
make_tabular(x, dd=2) {
  nr <- nrow(x)
  nc <- ncol(x)
  xf <- format(x, digits=dd)  # make them strings formatted as desired
  xf <- cbind(rownames(xf), xf)  # add the rownames as a column
  xf <- rbind(c(names(dimnames(xf))[1], dimnames(xf)[[2]]), xf) # add the colnames as the top row
  
  write.table(format(x, digits=2), quote=F, sep="  &  ", col.names=NA, eol="  \\\\  \n")
}
