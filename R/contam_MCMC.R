#' MCMC function for determing contamination probability and allele frequencies
#' 
#' Describtion
#' @param data 
#' @param inters 
#' @param ro_start
#' @param alpha
#' @param beta
#' @param lamda
#' 
#' @return Returns a list of two named components:
#' \describe{
#'  \item{prob_contam}{}
#'  \item{allele_freq}{}
#' }
#' @export
#' @examples
#'
#' @export
contam_MCMC<-function(data,inters,ro_start,alpha,beta,lamda){
  library(fullsniplings)
  snp_genos <- get_snp_genos(data) # Converts data from alleles of A,C,G,T to genotypes of 0, 1, or 2
  N <- ncol(snp_genos$mat) # Number of individuals
  L <- nrow(snp_genos$mat) # Number of loci
  ro <- rep(0,inters+1) # Creates array for ro values
  ro[1] <- ro_start 
  allele_f <- matrix(0,inters+1,L) # Matrix for allele frequencies
  allele_f[1,] <- rbeta(L,lamda,lamda)
  z <- matrix(0,inters,N) # Matrix for z's
  
  for(k in 1:inters){
  # update z
  prob <- full_z(snp_genos$mat,allele_f[k,],ro[k]) # full_z gives the prob of contamination for each individual
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
  for (i in 1:L){
    t_alpha <- 2*x2[i] + x1[i] + lamda # alpha parameter
    t_beta <- 2*x0[i] + x1[i] +lamda # beta paramenter
    allele_f[k+1,i] <- rbeta(1,t_alpha,t_beta) # allele frequency for 1 allele
  }

  # update ro
  sum_z <- sum(z[k,]) # total number of contaminated samples
  # updates ro with derived beta distribution
  p_alpha <- sum_z + alpha # alpha parameter
  p_beta <- N - sum_z + beta # beta parameter
  ro[k+1] <- rbeta(1,p_alpha,p_beta) # new ro value
  
  k = k + 1
  }
list(prob_contam = ro, allele_freq = allele_f)
}