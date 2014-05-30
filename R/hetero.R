#' Porportion of heterozygous loci
#' 
#' Computes the proportion of heterozygous loci of the sampled loci of one or 
#' more individuals.
#' @param gt A vector or matrix of genotypes.  They should be integers, 0,1, or 2 (or NA), where
#' the values record the number of copies f the "1" allele in the diploid geneotype. 
#' Note that if gt is a vector it can include the genotypes of more than one individual.  For example, 
#' there are L loci, the first L elements are loci of individual 1 and the next L
#' elements are loci of individual 2.  If gt is a matrix, each column represents an individual.
#' @param af A vector of allele frequencies for loci ordered as in gt giving the 
#' frequency of the 1 allele in the population.
#' @return Returns a list of one named component:
#' \describe{
#'  \item{p_heter}{Vector with number of elements equal to number of individuals sampled.
#'  Each element is proportion of heterozygous loci in one individuals sampled loci.}
#' }
#' @export
#' @examples
#' # run after random_gene to determine the heterozygousity of randomly generated genotypes
#' # determines proportion of heterozygous loci for 5 individuals
#' like <- likelihood(af = c(SNP1 = .1, SNP2 = .2))
#' genes <- random_gene(5,af)
#' hetero(genes$rclean,af)
hetero <- function(gt,af){
  nc <- length(af)
  gt <- matrix(gt, nrow=nc)
  heteros <- colSums (gt == 1)
  por <- heteros/nrow(gt)
  list(p_hetero = por)
}