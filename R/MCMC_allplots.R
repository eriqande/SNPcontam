#' @export
MCMC_allplots <- function(sim,rho,loci){
  MCMC_zplots(sim$z) # plots z values
  MCMC_alleleplot(sim$allele,rho,loci) # plots allele data
  MCMC_rhoplot(sim$rho) # plots rho data
  allele_table(sim$types,allele_df,sim$Lvals,sim$rhovals)
}