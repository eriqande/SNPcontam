

#' Chinook genotype data from 8031 individuals across much of the West Coast.
#' 
#' @docType data
#' @name swfsc_chinook_baseline
#' @usage data(swfsc_chinook_baseline)
#' @format This is a data frame with 8031 rows and 194 columns.  Each row contains
#' the genotype of a Chinook salmon included in the Southwest Fisheries Science
#' Center's Chinook baseline.  See Clemento et al.  The first four columns are 
#' descriptive of each fish:
#' \describe{
#'  \item{RepPop}{A catenation of the individual's reporting unit and population of origin.}
#'  \item{RepUnit}{The reporting unit of population from which the fish came.}
#'  \item{Pop}{The population from which the individual fish came.}
#'  \item{ID}{A unique identifier for every individual.  It turns out to be the RepPop catenated to a number.}
#' }
#' The remaining 190 columns are the alleles carried at each of the 95 loci.  This is in two-column
#' format, so every locus gets two adjacent columns.  The column headers have the names of the loci.  The
#' first occurrence of each locus is just the locus name and the second has a ".1" appended to it.
#' 
#' In general, alleles are coded as follows: A=1, C=2, G=3, T=4.  There are some loci that get an allele
#' code of 5, which refers to some data feature of some sort.  You can think of it as an allele.
#' There are only two alleles observed at each locus.  Missing data are denoted by NAs.  This data set
#' is described in detail in Clemento et al 2014 and is also available on Dryad (see references)
#' @keywords datasets
#' @references Clemento AJ, Crandall ED, Garza JC, Anderson EC (2014) Evaluation of a single nucleotide polymorphism baseline for genetic stock identification of Chinook Salmon (Oncorhynchus tshawytscha) in the California Current Large Marine Ecosystem. Fishery Bulletin 112(2-3): 112-130. doi:10.7755/FB.112.2-3.2
#' @references Clemento AJ, Crandall ED, Garza JC, Anderson EC (2014) Data from: Evaluation of a single nucleotide polymorphism baseline for genetic stock identification of Chinook Salmon (Oncorhynchus tshawytscha) in the California Current Large Marine Ecosystem. Dryad Digital Repository. doi:10.5061/dryad.574sv
NULL







#' Chinook genotype data from 1974 individuals sampled by CA Dept of Fish and Game at California Ports
#' 
#' @docType data
#' @name cdfg_port_samples
#' @usage data(cdfg_port_samples)
#' @format This is a data frame with 1974 rows and 190 columns.  Each row contains
#' the genotype of a Chinook salmon that was caught in the ocean and sampled at a port.
#' It was used for some analyses in Clemento et al.  The row name is the name given to the
#' fish.  The 190 columns are the alleles carried at each of the 95 loci.  This is in two-column
#' format, so every locus gets two adjacent columns.  The column headers have the names of the loci.  The
#' first occurrence of each locus is just the locus name and the second has a ".1" appended to it.
#' 
#' In general, alleles are coded as follows: A=1, C=2, G=3, T=4.  There are some loci that get an allele
#' code of 5, which refers to some data feature of some sort.  You can think of it as an allele.
#' There are only two alleles observed at each locus.  Missing data are denoted by NAs.
#' @keywords datasets
#' @references Clemento AJ, Crandall ED, Garza JC, Anderson EC (2014) Evaluation of a single nucleotide polymorphism baseline for genetic stock identification of Chinook Salmon (Oncorhynchus tshawytscha) in the California Current Large Marine Ecosystem. Fishery Bulletin 112(2-3): 112-130. doi:10.7755/FB.112.2-3.2
NULL
