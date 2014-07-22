#' MCMC function for determing contamination probability and allele frequencies
#' 
#' Description
#' @param MCMC A list that contains the outputs of the contam_MCMC function, including the rho output,
#' the allele frequency output, and the z value output from all sweeps.
#' @param burnin  The number of sweeps at the beginning of the MCMC that will be disregarded for analysis.
#' 
#' @return Returns a list of three named components:
#' \describe{
#'  \item{allele_pm}{A vector containing the posterior mean of the allele frequencies for each allele.}
#'  \item{rho_pm}{One value, which is the poterior mean of the contamination proportion.}
#'  \item{z_pm}{A vector containing the posterior mean of the z value for each individual.  A z_pm value
#'  1 at the ith index would indicate that the ith individual was indentified as contaminated in each sweep.}
#'}
#' @export
analyze_MCMC <- function(MCMC,burnin=100){
  alleles <- MCMC$allele_freq[-(1:burnin),]
  a_means <- colMeans(alleles)
  bottom <- apply(alleles, 2, quantile, probs = 0.05)
  top <- apply(alleles,2, quantile, probs = 0.95)
  rho_mean = mean(MCMC$prob_contam[-(1:burnin)])
  z <- MCMC$z[-(1:burnin),]
  z_pm <- colMeans(z)
  list(allele_pm = a_means, rho_pm = rho_mean, z_pm = z_pm)
}