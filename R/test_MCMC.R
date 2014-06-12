#' @export
test_MCMC <- function(sample_data, N, L, p, l, alpha, beta, lambda, inters, threshold){
  real_p <- p
  sim <- simulate_genos(N,L,p,l,sample_data)
  MCMC <- contam_MCMC(sim$geno,inters,rho_start=NULL,alpha,beta,lambda)
  
  #compare allele frequencies
  per_error <- (abs(colMeans(MCMC$allele_freq) - sim$afreqs)/sim$afreqs)*100
  error <- colMeans(MCMC$allele_freq) - sim$afreqs
  mean_error <- round(mean(abs(error)),3)
  plot(error)
  
  # calculate mean proportion contaminated
  p2 = round(mean(MCMC$prob_contam[-1]),3)
  
  # compare contaminated individuals
  contam_MCMC <- (colMeans(MCMC$z)>threshold)
  contam <- rep(0,N)
  contam[sim$contam_id] <- 1
  mistake <- contam_MCMC != contam
  false_n <- sum(mistake*contam)
  false_p <- sum(mistake) - false_n
  list(false_pos = false_p, false_neg = false_n, mean_error = mean_error, error_af = per_error, diff_af = error, mean_p = p2, MCMC_p = MCMC$prob_contam, MCMC_af = MCMC$allele_freq, afreqs = sim$afreqs, z = MCMC$z, contam = sim$contam_id)
  
}