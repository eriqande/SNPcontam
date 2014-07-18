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
# currently commented out because we don't use it.
#alle_tab <- allele_table(types = out_list_01$types, data = out_list_01$allele, Lvals = out_list_01$Lvals, rhovals = out_list_01$rhovals)


message("     making manuscript/tables/afreq_ci_coverage.tex       and   ")
message("            manuscript/tables/z_table_p_one_half.tex")

#### Produce a matrix for a the z-value table  ####
# we want a table that gives the fraction of truly contaminated samples with posterior prob 
# of contamination over 0.9 (or, more generally PPlim).  We also want the fraction of non-contaminated samples with that
# posterior prob >0.9 too.  We use the table function.
z <- out_list_01$z
PPlim <- 0.5

          
z_tab <- MCMC_ztable(z = out_list_01$z, PPlim = 0.5)



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

# here is a function to make output for tabular so we can customize more. 
make_tabular <- function(x, dd=2, leftcolhead, rightcolshead, outfile="", open_tabular = TRUE, close_tabular = TRUE, append = FALSE) {
  nr <- nrow(x)
  nc <- ncol(x)
  xf <- format(round(x, digits=dd), digits=dd)  # make them strings formatted as desired
  xf <- cbind(rownames(xf), xf)  # add the rownames as a column
  xf <- rbind(c(leftcolhead, colnames(x)), xf) # add the colnames as the top row
  
  
  # remove outfile if it exists
  if(append == FALSE) {if(file.exists(outfile)) file.remove(outfile)}
  
  # now make the enclosing tabular environment
  colstring <- paste(c("l", rep("r", ncol(x))), collapse="")
  if(open_tabular == TRUE) cat(paste("\\begin{tabular}{", colstring, "}", sep=""), "\n", file = outfile, append = T)
  
  # now, haggle with getting the multicolumn header on and written to file
  cat(paste("  &  \\multicolumn{", ncol(x), "}{c}{\\underline{", rightcolshead, "}} \\\\\n", sep=""), file = outfile, append = T)
  
  # now, the second line in xf should start with \hline
  xf[2,1] <- paste("\\hline", xf[2,1])
  write.table(format(xf, digits=2), quote=F, sep="  &  ", col.names=F, row.names=F, eol="  \\\\  \n", file = outfile, append = T)
  
  # finally, enclose the tabular environment if need be
  if(close_tabular == TRUE) cat("\\end{tabular} \n", file = outfile, append = T)
}



# Here we actually put the tables where they need to be to be inserted into the LaTeX document
make_tabular(a_table, leftcolhead = "$\\rho$~~~~~~~~~~", rightcolshead = "Number of Loci, $L$", outfile = "manuscript/tables/afreq_ci_coverage.tex")


# the z_table is a little lame at the moment.  I am just making it as two tables (a) and (b).  Hand customization would
# probably be better.  xtable doesn't allow imbedded multicolumns as far as I can tell, etc.  But maybe it is better
# this way since we need an extra digit for the false positive rates anyway.
ztab_file <- "manuscript/tables/z_table_p_one_half.tex"
cat("{\\bf (a)} Contaminated samples \n", file = ztab_file)
cat("\\begin{center}\n", file = ztab_file, append = T)
make_tabular(x = z_tab$TP[-1,], dd = 3, leftcolhead = "$\\rho~~~~~~~$", rightcolshead = "Number of Loci, $L$", outfile = ztab_file, append = T)
cat("\\end{center}\n", file = ztab_file, append = T)
cat("{\\bf (b)} Non-contaminated samples \n", file = ztab_file, append = T)
cat("\\begin{center}\n", file = ztab_file, append = T)
make_tabular(x = z_tab$FP, dd = 4, leftcolhead = "$\\rho~~~~~~~$", rightcolshead = "Number of Loci, $L$", outfile = ztab_file, append = T)
cat("\\end{center}\n", file = ztab_file, append = T)



