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
  
  # contamination population id
  which_contam <- which(z_pm >= 0.5)
  if(length(which_contam) >0){
  contam_pops <- lapply(which_contam, function(x) {if (z_pm[x] == 1) {MCMC$pop[-(1:(2*burnin)),x]} else{MCMC$pops[-(1:(2*burnin)),x][-c((2*which(MCMC$z[-(1:(burnin)),x] == 0)-1),(2*which(MCMC$z[-(1:(burnin)),x] == 0)))]}})
  # Contam population assignment
  opts <- combinations(P,2,1:P, repeats.allowed = TRUE) # all combinations of 2 populations
  t_opts <- nrow(opts) # number of combinations
  # Matrix of each individuals contamination population data. Two columns that are the two different combinations.
  mat_ind <- lapply(contam_pops, function(x) matrix(x,ncol=2, byrow = TRUE))
  # Number of times each of the combinations appears
  counts <- lapply(mat_ind, function(y) {sapply(1:t_opts, function(x) {sum(y[,1] == opts[x,1] & y[,2] == opts[x,2] | y[,1] == opts[x,2] & y[,2] == opts[x,1])})})
  # Pair with the maximum assignments for each individual
  max_pair <- sapply(counts, function(x) {a <- which(max(x) == x); if(length(a) > 1){a2 <- a[[1]]} 
                                            else{a2 <- a}; opts[a2,]})
  # Counts without the maximum
  counts2 <- lapply(counts, function(x) {a <- which(max(x) == x); if(length(a) > 1){a2 <- a[[1]]} 
                                           else{a2 <- a}; x[-a2]})
  # Pair with the second highest counts
  max_pair2 <- sapply(counts2, function(x) {a <- which(max(x) == x); if(length(a) > 1){a2 <- a[[1]]} 
                                              else{a2 <- a}; opts[a2,]})
  # Population Natmes
  populations = B_pops
  #Data frame wiht all of the name of the individual and the populations of origin according to the MCMC
  contam_df <- data.frame(id = which_contam, max_pair = matrix(populations[t(max_pair)],ncol=2), 
                            max_pair2 = matrix(populations[t(max_pair2)], ncol=2))
  } else{contam_df = NULL}
  
  list(z_pm = z_pm, rho_pm = rho_pm, pop_means = pop_means, pop_id = pop_id, contam_df = contam_df)
}
  
 