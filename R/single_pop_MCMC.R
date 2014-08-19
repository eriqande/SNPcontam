#' MCMC function for determing contamination probability and allele frequencies of single population sample
#' 
#' This function runs the single population Markov Chain Monte Carlo algorithm on a set SNP data to determine which samples
#' are contaminated, the proprotion of contaminated samples, and the allele frequencies of the SNP loci. The 
#' function:
#' \enumerate{
#' \item Changes the two column SNP data table into a L x N matrix where L is the number of loci, N is the 
#' number of individual samples, and genotypes are represented by 0s, 1s, and 2s.
#' \item Runs the MCMC algorithm on the data and stores each sweep.
#' \item Calculates the posterior means of the unknown values.
#' }
#' @param data A matrix containing the SNP data for all individuals.  The rows of the matrix are individual
#' samples and the columns are SNP data.  Each SNP loci has two alleles and therefore two columns.
#' @param locstart The index of the first column with genetic data.
#' @param inters  A number representing the total of number of interations for the MCMC model.
#' @param alpha The alpha parameter for the prior distribution of rho, the probability
#' of contamination.  The prior of rho is assumed to be a beta distribution.
#' @param beta The beta parameter for the prior distribution of rho, which is assumed to be
#' a beta distribution.
#' @param lambda The alpha and beta parameter for the prior distribution of the allele
#' frequencies at each locus, which is assumed to be a beta distribution.
#' @param burnin The number of sweeps at the beginning of the MCMC that will be disregarded for analysis.
#' 
#' @return Returns a list of two named components:
#' \itemize{
#'  \item{\strong{MCMC_data}: A list of three named components which comprise the raw output of the MCMC model:}
#'    \describe{
#'      \item{prob_contam}{A vector containing the rho value, which is the proportion of 
#'      contaminated samples for each interation. The vector is a length 1 plus the total number of iterations.}
#'      \item{allele_freq}{A matrix containing the allele frequencies at each locus for each
#'      iteration. The matrix has columns equal to the number loci and rows equal to the 1 plus
#'      total number of iterations.}
#'      \item{z}{A matrix containing the z value, which denotes contaminated status, for each individual
#'      and every iteration.  The matrix has columns equal to the number of individuals and rows equal to
#'      the total number of iterations.}
#'    }
#'  \item{\strong{analysis}: A list of three named components:}
#'      \describe{
#'        \item{allele_pm}{A vector containing the posterior mean of the allele frequencies for each allele.}
#'        \item{rho_pm}{One value, which is the poterior mean of the contamination proportion.}
#'        \item{z_pm}{A vector containing the posterior mean of the z value for each individual.  A z_pm value
#'        1 at the ith index would indicate that the ith individual was indentified as contaminated in each sweep.}}
#'      }
#' @export
#' @examples
#' data <- swfsc_chinook_baseline[1:500,]
#' # Run a short MCMC on the data
#' MCMC <- single_pop_MCMC(data, locstart = 5, inters = 100, burnin = 10)
#' names(MCMC)
#' names(MCMC$MCMC_data)
#' names(MCMC$analysis)
#' 
single_pop_MCMC <- function(data,locstart = 5, inters = 1000,alpha=0.5,beta=0.5,lambda=0.5, burnin = 100){
  change_to_ones_zeros_twos <- function(data, locstart){
    locs <- colnames(data)[locstart:ncol(data)]
  
    S <- data[, locs]  # grab just the genetic data
  
    uniq_alleles <- lapply(seq(1, ncol(S), 2), function(x) levels(factor(c(S[,x], S[, x+1]))))
    names(uniq_alleles) <- colnames(S)[seq(1, ncol(S), 2)]
  
    # drop loci that do not have two alleles
    Have2 <- sapply(uniq_alleles, function(x) length(x)==2)
    Have2rep <- rep(Have2, each=2)  # this can be used to extract the loci from M2 and B2 that have two alleles
  
    DropTheseLoci <- names(uniq_alleles)[!Have2]
    if(length(DropTheseLoci)>0) warning(paste("Dropping these loci that have > of < then 2 alleles:", paste(DropTheseLoci, collapse=", ")))
  
    S2 <- S[, Have2rep]
    uniq_alleles2 <- uniq_alleles[Have2]
  
    # now, make sure that the levels of S2's alleles are appropriately set
    for(i in seq(1, ncol(S2), 2)) {
      S2[, i] <- factor(S2[, i], levels = uniq_alleles2[[colnames(S2)[i]]])
      S2[, i+1] <- factor(S2[, i+1], levels = uniq_alleles2[[colnames(S2)[i]]])
    }
  
    tmp <- lapply(seq(1, ncol(S2), 2), function(i) {
      rowSums(
        cbind(
          as.integer(S2[, i]) - 1,   # the minus ones here make the alleles 0 and 1, rather than 1 and 2
          as.integer(S2[, i+1]) - 1
        )
      )}
    )
  
    prepared_data <- do.call(cbind, args = tmp)
    rownames(prepared_data) <- rownames(S2)
    colnames(prepared_data) <- colnames(S2)[seq(1, ncol(S2), 2)]
    prepared_data
  }
  
  prepared_data <- t(change_to_ones_zeros_twos(data = data,locstart = locstart))
  
  MCMC <- contam_MCMC(data = prepared_data,inters = inters,alpha = alpha, beta = beta,lambda = lambda)
  
  analyze <- analyze_MCMC(MCMC = MCMC, burnin = burnin)
  
  list(MCMC_data = MCMC, analysis = analyze)
  
} 