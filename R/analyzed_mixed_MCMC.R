#' @export
analyzed_mixed_MCMC <- function(MCMC, burnin = 100, B_pops){
  N <- ncol(MCMC$z)
  P <- length(B_pops)
  # Get posterior means for z and rho
  z_pm <- colMeans(MCMC$z[-(1:burnin),])
  rho_pm <- mean(MCMC$prob_contam[-(1:burnin)])
    
  # Posterior means for mixing proportions (pii)
  pii = colMeans(MCMC$mixing[-(1:burnin),])
    
  # Population assignment for clean populations
  clean_pops <- lapply(1:N, function(x) {MCMC$pops[-(1:(2*burnin)),x][2*which(MCMC$z[-(1:(2*burnin)),x] == 0) - 1]})
  tmp_pops <- lapply(clean_pops, function(x) {factor(x,level = 1:P)}) 
  freq <- lapply(tmp_pops, function(x) {as.numeric(table(x))})
  
  # fraction of assignments of each individual to each population
  pop_means <- do.call(what = rbind, args = lapply(freq, function(x) x))/nrow(MCMC$z)
  colnames(pop_means) <- B_pops
  
  # population id
  tmp_max <- matrix(sapply(1:N, function(x) max(pop_means[x,])),ncol=1)
  tmp_id <- sapply(1:N, function(x) colnames(pop_means)[which(max(pop_means[x,]) == pop_means[x,])])
  pop_id <- cbind(tmp_id,tmp_max)
  rownames(pop_id) <- rownames(pop_means)
  
  list(z_pm = z_pm, rho_pm = rho_pm, pop_means = pop_means, pop_id = pop_id)
}
  
 