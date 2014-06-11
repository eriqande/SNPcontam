#' Full Conditional Distribution of Contamination Indicator
#' 
#' Computes the probability of contamination for each individual sample from the full 
#' conditional distribution of the indicator value for contamination given the allele 
#' frequency of the 1 allele at each locus and the probability of contamination within 
#' the entire sample. 
#' @param genos A matrix containing the genotypes of all samples where 0s, 1s, and 2s 
#' are used to represent the three different genotype possiblities.  The number of
#' the number of loci, and the number of rows is equal to the number of individual samples.
#' @param theta A vector the ith element is allele frequency of the 1 allele
#' at the ith locus.
#' @param rho A number that is the probability that a sample is contaminated.
#' 
#' @return Returns a list of one named components:
#' \describe{
#'  \item{prob}{A vector, length equal to total number of samples, in which 
#'  the ith element is the probability that ith sample is contaminated.}
#' }
#' @export 
full_z <- function(genos,theta,rho){
  like <- likelihood(theta,genos) # use likelihood function to get probability of genotype given contamination and clean
  #consider taking log and then colSums so loop is avoided
  ln_clean <- log(like$clean) # take natural log of data
  ln_contam <- log(like$contam) # same for contaminated data
  clean <- (exp(colSums(ln_clean,na.rm=TRUE)))*(1 -rho) # probability of non contaminated sample
  contam <- (exp(colSums(ln_contam,na.rm=TRUE)))*rho # probability of contaminated samples
  # full conditional probability distribution of z indicator value of 1 (indicating contamination)
  # normalized by (clean + contam) so that total probability is 1
  p <- contam/(clean + contam) 
  return(p) # can just return value without list
}
  