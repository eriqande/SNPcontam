#' Likelihood ratio of genotype
#' 
#' Computes the likelihodd ratio of a genotype using the formula 
#' \eqn{\Lambda = -2*log(L(\theta_A|x)/L(\theta|x)$} where \eqn{L(\theta|x)} is the likelihood udner the assumption
#' that the sample is contaminated (DNA of two randomly sampled individuals is present)
#' and \eqn{L(\theta_0|x)} is the likelihood under the assumption that the sample is not
#' contaminated (DNA of only one individual is represented).
#' @param like_h A vector of likelihoods of the genotype at each loci undert he assumption.
#' that the sample is not contaminated. More than one individual can be included in lh.
#' For example, if there are L loci sampled, then the first L elements of lh are 
#' likelihoods of loci of individual 1 and the next L elements of lh are likelihoods 
#' of loci of individual 2.
#' @param like_hc Same as above but under the assumption that the sample is contaminated.
#' @param af A vector of the allele frequencies for loci giving the frequency of the 1
#' allele in a population.
#' @return Returns a list of one named component:
#' \describe{
#'  \item{rat}{Vector with number of elements equal to number of individuals sampled.
#'  Each element is the likelihood ratio of an individual's genotype.}
#' }
#' @export
#' @examples
#' # call it after likelihood function to likelihood ratios
#' # generates ratio for 5 individuals
#' like <- likelihood(af = c(SNP1 = .1, SNP2 = .2), c(0,0,0,1,1,1,1,1,1,2))
#' lratio(like$clean, like$contam,af)
lratio <- function(like_h, like_hc,af) {
  nc <- length(af) # nc is number of loci
  lh <- matrix(like_h, nrow=nc) # makes the likehoods into a matrix in case they are input as vectors
  lhc <- matrix(like_hc, nrow=nc) # does same as above for contaminated likelihoods
  lc <- apply(lhc,2,prod) # multplies all elements in the columns, corresponds with individual likelihood
  lnc <- apply(lh,2,prod) # does same as above but for contaminated individuals
  rat <- 2*log(lc/lnc) # calculates test statistic = likelihood ratio
  list(ratio = rat) # output ratio as list
}