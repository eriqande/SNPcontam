#' MCMC function for determing contamination porpotion and allele frequencies
#' 
#' This function runs an MCMC algorithm on the SNP genotype data of a group of samples to estimate the contamation
#' proportion and the allele frequencies and identify the contaminated samples.
#' @param data A L x N matrix containing the genotype of data of individuals in the form of 0s,1s, and 2s.
#' N is the number of individuals, and L is the number of loci.
#' @param inters  A number representing the total of number of interations for the MCMC model.
#' @param alpha The alpha parameter for the prior distribution of rho, the probability
#' of contamination.  The prior of rho is assumed to be a beta distribution.
#' @param beta The beta parameter for the prior distribution of rho, which is assumed to be
#' a beta distribution.
#' @param lambda The alpha and beta parameter for the prior distribution of the allele
#' frequencies at each locus, which is assumed to be a beta distribution.
#' 
#' @return Returns a list of three named components:
#' \describe{
#'  \item{prob_contam}{A vector containing the rho value, which is the proportion of 
#'  contaminated samples for each interation. The vector is a length 1 plus the total number of iterations.}
#'  \item{allele_freq}{A matrix containing the allele frequencies at each locus for each
#'  iteration. The matrix has columns equal to the number loci and rows equal to the 1 plus
#'  total number of iterations.}
#'  \item{z}{A matrix containing the z value, which denotes contaminated status, for each individual
#'  and every iteration.  The matrix has columns equal to the number of individuals and rows equal to
#'  the total number of iterations.}
#' }
#' @export
contam_MCMC<-function(data,inters = 1000,alpha=0.5,beta=0.5,lambda=0.5){
  N <- ncol(data) # Number of individuals
  L <- nrow(data) # Number of loci
  rho <- rep(0,inters+1) # Creates array for rho values
  rho[1] <- rbeta(1,alpha,beta)  
  
  allele_f <- matrix(0,inters+1,L) # Matrix for allele frequencies
  allele_f[1,] <- rbeta(L,lambda,lambda)
  z <- matrix(0,inters,N) # Matrix for z's
  
  for(k in 1:inters){
  # update z
  prob <- full_z(data,allele_f[k,],rho[k]) # full_z gives the prob of contamination for each individual
  z[k,] <- c(runif(N) <= prob)*1 # sets zi's to be 1 or 0 dependent on prob of contamination
  
  # update allele frequency
  # only use non-contaminated samples to calculate allele frequency
  # gene_x has 1s at indices where z is 0 and the genetype is x
  gene_0 <- (1 - matrix(rep(z[k,],L),nrow=L,byrow=TRUE))*(data==0)
  gene_1 <- (1 - matrix(rep(z[k,],L),nrow=L,byrow=TRUE))*(data==1)
  gene_2 <- (1 - matrix(rep(z[k,],L),nrow=L,byrow=TRUE))*(data==2)
  # sum of 1s of rows of gene_x gives total number of total number of genotype x at each locus
  x0 <- rowSums(gene_0,na.rm=TRUE) # total 0 genotype 
  x1 <- rowSums(gene_1,na.rm=TRUE) # total 1 genotype
  x2 <- rowSums(gene_2,na.rm=TRUE) # total 2 genotype
  # updates allele frequency using derived beta distribution
    t_alpha <- 2*x2 + x1 + lambda # alpha parameter
    t_beta <- 2*x0 + x1 + lambda # beta paramenter
    allele_f[k+1,] <- rbeta(L,t_alpha,t_beta) # allele frequency for 1 allele

  # update rho
  sum_z <- sum(z[k,]) # total number of contaminated samples
  # updates rho with derived beta distribution
  p_alpha <- sum_z + alpha # alpha parameter
  p_beta <- N - sum_z + beta # beta parameter
  rho[k+1] <- rbeta(1,p_alpha,p_beta) # new rho value

  }
list(prob_contam = rho, allele_freq = allele_f, z=z)
}