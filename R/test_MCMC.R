#' @export
test_MCMC <- function(sample_data, N, L, p, l, alpha=0.5, beta=0.5, lambda=0.5, inters, burnin=100){
  real_p <- p
  sim <- simulate_genos(N,L,p,l,sample_data)
  MCMC <- contam_MCMC(sim$geno,inters,rho_start=NULL,alpha,beta,lambda)
  
  #compare allele frequencies
  alleles <- MCMC$allele_freq[-(1:burnin),]
  a_means <- colMeans(alleles) 
  bottom <- apply(alleles, 2, quantile, probs = 0.05)
  top <- apply(alleles,2, quantile, probs = 0.95)
  
  # calculate mean proportion contaminated
  p2 = mean(MCMC$prob_contam[-(1:burnin)])
  
  # compare contaminated individuals
  z <- MCMC$z[-(1:burnin),]
  z_pm <- colMeans(MCMC$z)
  z_id <- rep(0,N)
  z_id[sim$contam_id] = 1
  
  list(afreqs = sim$afreqs, bottom_int = bottom, top_int = top , alle_pm = a_means ,z_pm = z_pm , z_id = z_id, rho_pm = p2 )
  
}