# this is a script to perform the simulations and analyses, make all the figures
# then typset all the documents for this paper.  


#### Load Packages and Check that We Are in the Right Directory ####
library(knitr)

# this script should be run from the directory that includes the directories "simulations", "supplements", 
# and "manuscript"
if(!all(file.exists("simulations", "manuscript", "supplements")))  {
  stop("You must run 01_simulation_1.R in directory that includes: \"simulations\", \"manuscript\", \"supplements\"")
}


#### Do Simulation Analyses  ####
# The following line is commented out because we have the output: out_list_01.rda, already.
# If you wanted to regenerate that output, uncomment it, but plan on using a lot of cores
# to get the job done quickly.
#source("simulations/01_simulation_1.R")  


#### Create the Figures and Tables  ####
# note that these scripts place the figures and tables where they
# need to be relative to the LaTeX documents.
source("simulations/02_makecharts_1.R")
source("simulations/03_maketables_1.R")


#### Crop Figures As Needed  ####
# before we start compiling up the documents, we are going to 
# run pdfcrop on all of our figures
# this is a quick funtion to run pdfcrop on all the pdf files in a directory
pdfcrop_dir <- function(path = "manuscript/images/", pattern = "*.pdf", exclude = "banner.pdf") {
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
curdir <- getwd()
setwd("manuscript/")
knit2pdf("snp-contam-main.Rnw")
setwd(curdir)


#### Move Manuscripts and Supplements to an Output Directory ####
# not implemented yet.
