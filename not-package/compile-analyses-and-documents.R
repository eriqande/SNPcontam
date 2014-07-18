# this is a script to perform the simulations and analyses, make all the figures
# then typset all the documents for this paper.  

# NOTE: By default this will not re-do the lengthy simulations.  Rather it will
# just read in stored results from the files ./simulations/out_list_01.rda and
# ./simulations/out_list_02.rda.   You can change variables below to force the 
# simulations to be re-done, alter that, but be warned that it will try to run things
# on 20 processors with mclapply and will take a long time...
REDO_SIMULATION_1 <- FALSE
REDO_SIMULATION_2 <- FALSE


#### Check to see if necessary packages are installed, and install if necessary ####
message("Checking that required packages are installed")
packages <- c("ggplot2", "grid", "gridExtra", "gtools", "knitr", "parallel", "reshape2", "xtable")
need_these <- setdiff(packages, rownames(installed.packages()))

if(length(need_these) > 0) {
   message("We are going to attempt to install the following packages (which are not yet installed):")
   message(need_these)
   install.packages(need_these)
}



#### Load Packages and Check that We Are in the Right Directory ####
library(knitr)

# this script should be run from the directory that includes the directories "simulations", "supplements", 
# and "manuscript"
if(!all(file.exists("simulations", "manuscript", "supplements")))  {
  stop("You must run compile-analyses-and-documents.R in directory that includes: \"simulations\", \"manuscript\", \"supplements\"")
}




#### Do Simulation And Make Figures for Single Source Population Model  ####
# The following line is commented out because we have the output: out_list_01.rda, already.
# If you wanted to regenerate that output, uncomment it, but plan on using a lot of cores
# to get the job done quickly.

if(REDO_SIMULATION_1 == TRUE) {
  message("Starting to re-do simulation 1.  This could take a while!")
  source("simulations/01_simulation_1.R")    # this creates the file ./simulations/out_list_01.rda 
}

# these scripts place the figures and tables where they need to be relative to the LaTeX documents.
message("Making figures from results from simulation 1")
source("simulations/02_makecharts_1.R")

message("Making tables from results from simulation 1")
source("simulations/03_maketables_1.R")




#### Do Simulation And Make Figures for Mixed Fishery Model  ####
# The following line is commented out because we have the output: out_list_01.rda, already.
# If you wanted to regenerate that output, uncomment it, but plan on using a lot of cores
# to get the job done quickly.

if(REDO_SIMULATION_2 == TRUE) {
  message("Starting to re-do simulation 2.  This could take a while!")
  set.seed(10)  # for reproducibility of results
  source("simulations/04_simulation_2.R")    # this creates the file ./simulations/out_list_02.rda 
}

# this script places the figures and tables where they need to be relative to the LaTeX documents.
message("Making figures and tables from the results of simulation 2")
source("simulations/05_makecharts_tables_2.R")





#### Crop Figures As Needed  ####
message("Using pdfcrop to crop PDF figures.  This requires a decent LaTeX installation that has pdfcrop")
# before we start compiling up the documents, we are going to 
# run pdfcrop on all of our figures
# this is a quick funtion to run pdfcrop on all the pdf files in a directory
pdfcrop_dir <- function(path = "manuscript/images/", pattern = "*.pdf", exclude = c("banner.pdf", "DAG_1.pdf", "DAG_2.pdf")) {
  img <- dir(path = path, pattern = pattern)
  img <- img[!(img %in% exclude)]
  for(i in img) {
    system(paste("pdfcrop", file.path(path, i), file.path(path, i)))
  }
}
# now, do the cropping:
pdfcrop_dir()
pdfcrop_dir(path = "supplements/images", exclude = NULL)




#### Compile the Manuscript  ####
message("Knitting the manuscript to PDF.  This requires a decent LaTeX installation")
curdir <- getwd()
setwd("manuscript/")
knit2pdf("snp-contam-main.Rnw")
setwd(curdir)


#### Move Manuscripts and Supplements to an Output Directory ####
# not implemented yet.
