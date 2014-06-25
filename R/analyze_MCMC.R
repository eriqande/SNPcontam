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