#' MCMC function for determing contamination proportions and mixing proportions of a mixture sample
#' 
#' Description
#' @param data A L x N matrix containing the genotype of data of individuals in the form of 0s,1s, and 2s.
#' N is the number of individuals, and L is the number of loci.
#' @param contam_data A P*P x N matrix containing the likelihood that each individual orginates from each
#' combination of the 2 populations from P, assuming the individual is contaminated.  P is the number of
#' populations and is the number of individuals.
#' @param clean_data A P x N matrix containing the likelihood that each individual originates from each
#' of P populations, assuming the individual is not contaminated.
#' @param alpha The alpha parameter for the prior distribution of rho, the probability
#' of contamination.  The prior of rho is assumed to be a beta distribution.
#' @param beta The beta parameter for the prior distribution of rho, which is assumed to be
#' a beta distribution.
#' @param inters Total number of sweeps used in the MCMC model.
#' 
#' @return Returns a list of four named components:
#' \describe{
#'  \item{prob_contam}{A vector containing the rho value, which is the proportion of 
#'  contaminated samples for each interation. The vector is a length of 1 plus the total number of iterations.}
#'  \item{pops}{A matrix the u, which denote the population or populations of origin, for each sweep.
#'  u is a 2*inters x N matrix, and each pair of rows represents the population identification for each sweep.
#'  If an individual has u values of (any P, 0), than it is a non contaminated individual and the first value is
#'  its population of origin.  If an individual has two u values, it is contaminated and the two values 
#'  represent the source of contamination.}
#'  \item{z}{A matrix containing the z value, which denotes contaminated status, for each individual
#'  and every iteration.  The matrix has columns equal to the number of individuals and rows equal to
#'  the total number of iterations.}
#'  \item{mixing}{A matrix containing the mixing proportions of each population for each sweep.  The matrix
#'  has number of rows equal to the total interations plus one and columns equal to the number of populations.}
#' }
#' @export
#' @examples
#' 
#' # Creates data for the mixed_MCMC function
#' zeros <- t(matrix(c(40,160,50,80,100,20,30,15,140),nrow=3))
#' ones <- 200 - zeros
#' genos <- t(matrix(c(1,2,0,2,1,0,1,0,0,1,2,2,1,1,1),nrow=5))
#' data <- matrix(c(1,2,1,0,2,0,2,0,0,0,1,2,0,0,1),nrow=3)
#' clean_data <- P_likelihood(zeros,ones,genos,.5)
#' contam_data <- Pcontam(zeros,ones,genos,.5) 
#' # Runs the MCMC on the data for 5 interations
#' mixed_MCMC(data, contam_data, clean_data, inters = 5)
mixed_MCMC <- function(data, contam_data, clean_data, alpha=.5, beta=.5, inters){
  
  cbs <- nrow(contam_data) # Number of combinations of populations
  P <- nrow(clean_data) # Number of populations
  N <- ncol(data) # Number of individuals
  L <- nrow(data) # Number of loci
  
  # Sets up rho matrix and gets first rho value
  rho <- rep(0,inters+1)
  rho[1] <- rbeta(1,alpha,beta)
  
  # Sets up u matrix
  u <- matrix(0,N,nrow=inters*2)
  
  # Sets up pi matrix and gets first pi values
  pii <- matrix(0,inters+1,P)
  gm0 <- rgamma(P,shape=rep(1/P,P))
  pii[1,] <- gm0/sum(gm0)
  
  # Sets up z matrix
  z <- matrix(0,inters,N)
  
  # Creates list containing all combinations of populations
  combos <- matrix(0,2,P*P)
  pops <- 1:P
  k <- 0
  for (p1 in pops) for (p2 in pops){
    k <- k + 1
    combos[ ,k] <- c(p1,p2)
  }
  
  for(i in 1:inters){
    # propabilities of an individual coming from combinations of populations given the mixture proportions
    c_pi <- rep(pii[i,], each=P)*rep(pii[i,], times=P)
  # update z
    probs_c <- contam_data*c_pi # prob of originating from combinations given genotypes and mixture proportions
    probs0 <- clean_data*pii[i,] # prob of originating from each singular population given genotypes and mixture proportions
    pcontam <- colSums(probs_c*rho[i]) #prob of individuals being contaminated
    pclean <- colSums(probs0*(1-rho[i])) # prob of individuals being noncontamined
    prob <- pcontam/(pclean + pcontam) # normalize probability so it adds to one
    z[i,] <- c(runif(N) <= prob)*1 # chooses new z based on prob of being contaminated
  
  # update u
    z_c <- which(z[i,] %in% 1) # indices of contaminated individuals
    z0 <- which(z[i,] %in% 0) # indices of clean individuals
    #contaminated
    if (length(z_c) != 0){
      if (length(z_c) == 1){
        uprobs_c <- probs_c[,z_c] # prob of originating from each combination
        pops_c <- sample((1:cbs), size=1, prob = uprobs_c) # samples from combinations based on uprobs_c
        }
      else {
        uprobs_c <- apply(probs_c[,z_c],2,list) # prob of orginating from each combo for each individual in z_c
        pops_c <- sapply(uprobs_c, function(x) {sample((1:cbs),size=1, prob = unlist(x))}) # samples for each individual in z_c
        }
      u[(2*i-1):(2*i), z_c] <- combos[pops_c] # updates u for contaminated individuals
    }
    #clean
    if (length(z0) != 0){
      if (length(z0) == 1){
        uprobs0 = probs0[,z0] # prob of orginating from each population
        pops0 <- sample((1:P),size=1, prob = uprobs0) # samples from all populations given uprobs0
      }
      else {
        uprobs0 <- apply(probs0[,z0],2,list) # prob of orginating from each population for each individual in z0
        pops0 <- sapply(uprobs0, function(x) {sample((1:P), size=1, prob = unlist(x))}) # samples population for each individual in z0
      }
      u[(2*i-1),z0] <- pops0 # updates u for clean individuals
    }  
  
  # update pi
    clean_us <- u[1,][u[2,]==0] # the population identification for only clean individuals
    utmp <- factor(clean_us, level=(1:P)) # step before using table so table will include values that have frequency of 0
    freqs <- as.numeric(table(utmp)) # gets frequencies of each population 
    xi <- (freqs + 1/P) # parameters for the dirichlet distribution for pi
    # samples from dirichlet distribution using gamma distribution and updats pii
    gm <- rgamma(P,shape=(xi)) 
    pii[i+1,] <- gm/sum(gm)
  
  # update rho
  sum_z <- sum(z[i,]) # total number of contaminated samples
  # updates rho with derived beta distribution
  p_alpha <- sum_z + alpha # alpha parameter
  p_beta <- N - sum_z + beta # beta parameter
  rho[i+1] <- rbeta(1,p_alpha,p_beta)
  }
list(prob_contam = rho, pops = u, z = z, mixing = pii)
}