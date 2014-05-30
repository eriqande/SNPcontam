#' Random generation of genotypes given genotype probabilities
#' 
#' Generates genotypes (using 0, 1, and 2) given the probabilities of 
#' the genotype at each locus given that the sample is non-contaminated
#' (DNA of only one individual is represented) and given that the sample is
#' contaminated (DNA of two randomly sampled individuals is present).
#' @param n An integer representing the desired number of genotypes to be generated
#' @param af A vector of the allele frequencies for loci giving the frequency of the 1
#' allele in a population.
#' @return Returns a list of two named components:
#' \describe{
#'  \item{rclean}{Matrix with number rows equal to the number of loci (
#'  the length of af) and n columns representing n simulated genotypes.
#'  The given genotypes represent non-contaminated samples.}
#'  \item{rcontam}{Same as above but genotypes represent contamination samples.}
#' }
#' @export
#' @examples
#' # generates 5 different genotypes
#' random_gene(5,af = c(SNP1 = .1, SNP2 = .2))
random_gene <- function(n,af) {
  L <- likelihood(af) # calls likelihood function which contans genotype probability
  l = length(af) # number of loci
  # for g matrices: columns are individuals and rows are loci
  g1 <- matrix(0, l, n) # matrix for non contamined
  g2 <- g1 # matrix for contamined individuals
  # loop cycles through each locus
  for (i in 1:l){
    # stores genotype at locus in all individuals
    #L$clean_prob is matrix fo genotype probabilities or each locus
    g1[i,] = sample(0:2, n, TRUE, prob = L$clean_prob[,i])
    # same as above but for contaminated
    g2[i,] = sample(0:2, n, TRUE, prob = L$contam_prob[,i]) 
  }
  list(rclean = g1, rcontam = g2) # output the two random groups of genotypes as list
}