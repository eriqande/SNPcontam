#' Probability of genotype given contamination and noncontamination
#' 
#' Computes the probability of a genotypes (0, 1, or 2) given the 
#' allele frequencies in the population and the assumption that it is 
#' non-contaminated (DNA of only one individual is represented) or is 
#' contaminated (DNA of two randomly sampled individuals is present)
#' @param af A vector of allele frequencies for loci ordered as in gt
#' giving the freqency of the 1 allele in the population.  Note that gt
#' could be much longer than af. For example, if gt is a vector where say the first 
#' L elements are the L loci at individual 1 and the next L are the loci at 
#' individual 2 and so on.  In that case, af will recycle as it should.
#' @param gt Vector of genotypes.  They should be integers, 0, 1, or 2
#' (or NA), where the values record the number of copies of the "1" allele
#' in the diploid genotype.  This can be omitted, in which case the just
#' the genotype probability matrices are returned
#' @return Returns a list of four named components:
#' \describe{
#'  \item{contam_prob}{Matrix with three rows (corresponding the genotypes 0, 1, and 2)
#'  and number of columns equal to length of af.  Given probabilities of the three possible
#'  genotypes for each SNP under the hypothesis of contamination.}
#'  \item{clean_prob}{Same as above but for no contamination.}
#'  \item{contam}{Vector of genotype probabilities under the assumption of contamination}
#'  \item{clean}{Vector of genotype probabilities under the assumption of no contamination}
#' }
#' @export
#' @examples
#' # call it with just allele frequencies to get the genotype probability for two SNPs.
#' # note that af can have names
#' likelihood(af = c(SNP1 = .1, SNP2 = .2))
#' 
#' # call it with af and also genotypes you want to evaluate probabilities for
#' likelihood(af = c(SNP1 = .1, SNP2 = .2), c(0,0,0,1,1,1,1,1,1,2))
#' 
#' # for illustration, get the contam probs and turn them into a matrix
#' # of values where the individuals are columns and the rows are SNPs
#' matrix(likelihood(af = c(SNP1 = .1, SNP2 = .2), c(0,0,0,1,1,1,1,1,1,2))$contam, nrow = 2)
likelihood <- function(af, gt = integer(0)) {
  
  # make the list to return and set values to NULL
  ret <- list(contam_prob = NULL, clean_prob = NULL, contam = NULL, clean = NULL)
  
  #calculate probablitiy for noncontaminated
  g0 <- (1-af)^2; g1 <- 2*af*(1-af); g2 <- af^2
  ret$clean_prob <- rbind(g0,g1,g2)
  colnames(ret$clean_prob) <- names(af)
  
  #calculate probability for contaminated
  gc0 <- (1-af)^4; gc1 <- 1 - (af^4 + (1-af)^4); gc2 <- af^4
  ret$contam_prob <- rbind(gc0, gc1, gc2)
  colnames(ret$contam_prob) <- names(af)
  
  #make vector of likelihood for each non contaminated gene, if there is a gt vector of length >0.
  if(length(gt)>0){
    a <- gt==0; b <- gt==1; c <- gt==2
    ret$clean <- g0*a + g1*b + g2*c
    
    #make vector of likelihood for each contaminated gene
    ac <- gt==0; bc <- gt==1; cc <- gt==2
    ret$contam <- gc0*ac + gc1*bc + gc2*cc
  }  
  
  ret
}
