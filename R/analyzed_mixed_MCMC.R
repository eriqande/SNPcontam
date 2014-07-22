#' Analyzed the output of an MCMC algorithm for a mixture model.
#' 
#' This function takes the raw output of the MCMC algorithm function \code{mixed_MCMC} and calulates the 
#' posterior means of the unknown values and performs the population assignment.
#' @param MCMC A list that contains the outputs of the contam_MCMC function, including the rho output,
#' the population assignment output, the mixture proportion output, and the z value output from all sweeps.
#' @param burnin The number of sweeps at the beginning of the MCMC that will be disregarded for analysis.
#' @param B_pops A factor vector giving the population of origin of each individual in the Baseline.  
#' This should have levels which are ordered the way they should be.
#' 
#' @return Returns a list of four named components:
#' \describe{
#'  \item{z_pm}{A vector containing the posterior mean of the z value for each individual.  A z_pm value
#'  1 at the ith index would indicate that the ith individual was indentified as contaminated in each sweep.}
#'  \item{rho_pm}{One value, which is the poterior mean of the contamination proportion.}
#'  \item{pop_means}{A N x P matrix where N is the number of individual samples and P is the number of populations
#'  in the baseline. Each element is the fraction of times that a certain individual was assigned to a specific
#'  population.}
#'  \item{pop_id}{A N x 2 matrix that contatins the populaton of assignment for each individual and the fraction
#'  of sweeps in which the individual was assigned to this population.  The rownames of the matrix are the names of
#'  each individual.}
#' }
#' @export
#' @examples
#' # grab on MCMC run on simulated data
#' load("~/Hollings/SNPcontam/not-package/MCMC.rda")
#' test_MCMC <- MCMC[[1]][[1]]$output$MCMC
#' analysis <- analyzed_mixed_MCMC(test_MCMC, B_pops = levels(swfsc_chinook_baseline$RepPop))
#' names(analysis)
#' 
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
  
 