#' MCMC function for determing contamination probability and allele frequencies
#' 
#' Description
#' @param data A matrix containing the genotype of data of individuals in the form of 
#' A,C,G, and T.
#' @param inters  A number representing the total of number of interations for the MCMC model.
#' @param rho_start The start value for the probability of contamination within the entire
#' sample.
#' @param alpha The alpha parameter for the prior distribution of rho, the probability
#' of contamination.  The prior of rho is assumed to be a beta distribution.
#' @param beta The beta parameter for the prior distribution of rho, which is assumed to be
#' a beta distribution.
#' @param lambda The alpha and beta parameter for the prior distribution of the allele
#' frequencies at each locus, which is assumed to be a beta distribution.
#' 
#' @return Returns a list of two named components:
#' \describe{
#'  \item{prob_contam}{A vector containing the rho value, which is the probability of 
#'  contamination, for each interation. The vector has 1 plus the total number of iterations.}
#'  \item{allele_freq}{A matrix containing the allele frequencies at each locus for each
#'  iteration. The matrix has columns equal to the number loci and rows equal to the 1 plus
#'  total number of iterations.}
#' }
#' @export
contam_MCMC<-function(data,inters,rho_start,alpha,beta,lambda){
  library(fullsniplings)
  snp_genos <- get_snp_genos(data) # Converts data from alleles of A,C,G,T to genotypes of 0, 1, or 2
  N <- ncol(snp_genos$mat) # Number of individuals
  L <- nrow(snp_genos$mat) # Number of loci
  rho <- rep(0,inters+1) # Creates array for rho values
  rho[1] <- rho_start 
  allele_f <- matrix(0,inters+1,L) # Matrix for allele frequencies
  allele_f[1,] <- rbeta(L,lambda,lambda)
  z <- matrix(0,inters,N) # Matrix for z's
  
  for(k in 1:inters){
  # update z
  prob <- full_z(snp_genos$mat,allele_f[k,],rho[k]) # full_z gives the prob of contamination for each individual
  z[k,] <- c(runif(N) <= prob$prob)*1 # sets zi's to be 1 or 0 dependent on prob of contamination
  
  # update allele frequency
  # only use non-contaminated samples to calculate allele frequency
  # gene_x has 1s at indices where z is 0 and the genetype is x
  gene_0 <- (1 - matrix(rep(z[k,],L),nrow=L,byrow=TRUE))*(snp_genos$mat==0)
  gene_1 <- (1 - matrix(rep(z[k,],L),nrow=L,byrow=TRUE))*(snp_genos$mat==1)
  gene_2 <- (1 - matrix(rep(z[k,],L),nrow=L,byrow=TRUE))*(snp_genos$mat==2)
  # sum of 1s of rows of gene_x gives total number of total number of genotype x at each locus
  x0 <- rowSums(gene_0,na.rm=TRUE) # total 0 genotype 
  x1 <- rowSums(gene_1,na.rm=TRUE) # total 1 genotype
  x2 <- rowSums(gene_2,na.rm=TRUE) # total 2 genotype
  # updates allele frequency using derived beta distribution
    t_alpha <- 2*x2 + x1 + lambda # alpha parameter
    t_beta <- 2*x0 + x1 + lambda # beta paramenter
    allele_f[k+1,] <- rbeta(1,t_alpha,t_beta) # allele frequency for 1 allele

  # update rho
  sum_z <- sum(z[k,]) # total number of contaminated samples
  # updates rho with derived beta distribution
  p_alpha <- sum_z + alpha # alpha parameter
  p_beta <- N - sum_z + beta # beta parameter
  rho[k+1] <- rbeta(1,p_alpha,p_beta) # new rho value
  
  k = k + 1
  }
list(prob_contam = rho, allele_freq = allele_f, z=z)
}