#' MCMC function for determing contamination proportion and contaminated samples of mixture sample
#' 
#' This function runs the multiple population Markov Chain Monte Carlo algorithm on a set SNP data to determine which samples
#' are contaminated, the proprotion of contaminated samples, and the mixing proportions of the different populations.
#' \enumerate{
#' \item Changes the two column SNP data table into a L x N matrix where L is the number of loci, N is the 
#' number of individual samples, and genotypes are represented by 0s, 1s, and 2s.
#' \item Runs the MCMC algorithm on the data and stores each sweep.
#' \item Calculates the posterior means of the unknown values and records the population identification.
#' }
#' @param baseline A matrix containing the SNP data for all individuals in the baseline.The rows of the matrix 
#' are individual samples and the columns are the SNP data,  Each SNP loci has two alleles and therefore 
#' two columns.
#' @param mixture A matrix containing the SNP data for all individuals to be tested for contamination.  
#' The matrix is in the same format as the baseline data.
#' @param B_locstart The index of the first column with genetic data in the baseline.
#' @param B_pops A factor vector giving the population of origin of each individual in the Baseline.  This should have levels 
#' which are ordered the way they should be.
#' @param M_locstart The index of the first column with genetic data in the mixture data.
#' @param alpha The alpha parameter for the prior distribution of rho, the probability
#' of contamination.  The prior of rho is assumed to be a beta distribution.
#' @param beta The beta parameter for the prior distribution of rho, which is assumed to be
#' a beta distribution.
#' @param lambda The alpha and beta parameter for the prior distribution of the allele
#' frequencies at each locus, which is assumed to be a beta distribution.
#' @param inters  A number representing the total of number of interations for the MCMC model.
#' @param burnin The number of sweeps at the beginning of the MCMC that will be disregarded for analysis.
#' @param contamination Logical variable.  If contamination = FALSE, then the MCMC model will only update 
#' the mixing proportions and the population identity, assuming no contamination is present.
#' 
#' @return Returns a list of two named components:
#' \itemize{
#'  \item{\strong{MCMC_data}: A list of four named components which comprise the raw output of the MCMC model:}
#'    \describe{
#'      \item{prob_contam}{A vector containing the rho value, which is the proportion of 
#'      contaminated samples for each interation. The vector is a length of 1 plus the total number of iterations.}
#'      \item{pops}{A matrix the u, which denote the population or populations of origin, for each sweep.
#'      If contamintion=TRUE, u is a 2*inters x N matrix, and each pair of rows represents the population identification for each sweep.
#'      If an individual has u values of (any P, 0), than it is a non contaminated individual and the first value is
#'      its population of origin.  If an individual has two u values, it is contaminated and the two values 
#'      represent the source of contamination.  If contamination=FALSE, u is a intersxN matrix, and each row represents the popultion
#'      identification for each sweep.}
#'      \item{z}{A matrix containing the z value, which denotes contaminated status, for each individual
#'      and every iteration.  The matrix has columns equal to the number of individuals and rows equal to
#'      the total number of iterations.}
#'      \item{mixing}{A matrix containing the mixing proportions of each population for each sweep.  The matrix
#'      has number of rows equal to the total interations plus one and columns equal to the number of populations.}}

#'    }
#'  \item{\strong{analysis}: A list of three named components:}
#'      \describe{
#'        \item{z_pm}{A vector containing the posterior mean of the z value for each individual.  A z_pm value
#'        1 at the ith index would indicate that the ith individual was indentified as contaminated in each sweep.}
#'        \item{rho_pm}{One value, which is the poterior mean of the contamination proportion.}
#'        \item{pop_means}{A N x P matrix where N is the number of individual samples and P is the number of populations
#'        in the baseline. Each element is the fraction of times that a certain individual was assigned to a specific
#'        population.}
#'        \item{pop_id}{A N x 2 matrix that contatins the populaton of assignment for each individual and the fraction
#'        of sweeps in which the individual was assigned to this population.  The rownames of the matrix are the names of
#'        each individual.}
#'      }
#' @export
#' @examples
#' mix <- sample(1:nrow(swfsc_chinook_baseline),100)
#' mixture <- swfsc_chinook_baseline[mix,]
#' baseline <- swfsc_chinook_baseline[-(mix),]
#' # Run a short MCMC on the data
#' MCMC <- multi_pop_MCMC(baseline=baseline,mixture=mixture, B_pops=baseline$RepPop, inters = 100, burnin = 10)
#' names(MCMC)
#' names(MCMC$MCMC_data)
#' names(MCMC$analysis)
 
multi_pop_MCMC <- function(baseline, mixture, B_locstart=5, B_pops, M_locstart=5, alpha=.5, beta=.5, lambda=.5, inters=1000, burnin=100,contamination=TRUE){
  require(gtools)
  prepared_data <- prepare_base_and_mix_for_mixed_MCMC(B = baseline, 
                                                       B_locstart = B_locstart, 
                                                       B_pops = B_pops, 
                                                       M = mixture, 
                                                       M_locstart)
  like_contam <- Pcontam(snp_zeroes = prepared_data$zeros,
                         snp_ones = prepared_data$ones, 
                         genos = prepared_data$mixmat,
                         lambda = lambda)
  like_clean <- P_likelihood(snp_zeroes = prepared_data$zeros,
                             snp_ones = prepared_data$ones,
                             genos = prepared_data$mixmat,
                             lambda = lambda)
  MCMC_raw_data <- mixed_MCMC(data = prepared_data$mixmat,
                              contam_data = like_contam, 
                              clean_data = like_clean,
                              alpha = alpha,
                              beta = beta, 
                              contamination = contamination,
                              inters = inters)
  MCMC_results <- analyzed_mixed_MCMC(MCMC = MCMC_raw_data, burnin = burnin, B_pops = levels(B_pops))
  list(MCMC_data = MCMC_raw_data, analysis = MCMC_results)
}