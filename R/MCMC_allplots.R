MCMC_allplots <- function(sim,rho,loci){
  MCMC_zplots(sim$z) # plots z values
  MCMC_alleleplot(sim$allele,rho,loci) # plots allele data
  MCMC_rhoplot(sim$rho, sim$rhovals) # plots rho data
  allele_table(sim$types,sim$allele,sim$Lvals,sim$rhovals)
}