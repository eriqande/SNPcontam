#' Porportion of heterozygous loci
#' 
#' Computes the proportion of heterozygous loci out of the total sampled loci for 
#' each individuals.
#' @param gt A vector or matrix of genotypes.  They should be integers, 0,1, or 2 
#' (or NA), where the values record the number of copies of the "1" allele in the 
#' diploid geneotype.  Note that if gt is a vector it can include the genotypes of more 
#' than one individual.  For example, there are L loci, the first L elements are loci of 
#' individual 1 and the next L elements are loci of individual 2.  If gt is a matrix, 
#' each column represents an individual.
#' @param af A vector of allele frequencies for loci ordered as in gt giving the 
#' frequency of the 1 allele in the population.
#' @return Returns a list of one named component:
#' \describe{
#'  \item{p_heter}{Vector with number of elements equal to number of individuals sampled.
#'  Each element is proportion of heterozygous loci in one individual's sampled loci.}
#' }
#' @export
#' @examples
#' # running the function on the genotype below returns that there are 4 heterozygotes
#' genotype <- c(1,1,0,0,1,2,1)
#' hetero(genotype,af = c(.5,.5,.2,.3,.6,.2,.1))
#' 
#' # run after random_gene to determine the heterozygousity of randomly generated genotypes
#' # determines proportion of heterozygous loci for 5 individuals
#' af = c(SNP1 = .1, SNP2 = .2)
#' like <- likelihood(af)
#' genes <- random_gene(5,af)
#' hetero(genes$rclean,af)
hetero <- function(gt,af){
  nc <- length(af) # number of loci
  # makes the genotypes into matrix with individuals in columns and genotype probs in case gt is vector
  geno <- matrix(gt, nrow=nc) 
  # vector containing number of heterozygous loci in each column/each individual
  heteros <- colSums (geno == 1)
  # calculates proportion of heterozygous loci for each individual
  por <- heteros/nrow(geno) 
  list(p_hetero = por) # output proportion as list
}