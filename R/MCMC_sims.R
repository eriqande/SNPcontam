#' @export
MCMC_sims <- function(sample_data, N, Lvals, rhovals, l=1, alpha=0.5, beta=0.5, lambda=0.5, inters=1000,n){
  # create list with all simulations
  library(ggplot2)
  library(grid)
  library(gridExtra)
  types <- list() 
  k <- 0
  num_Ls <- length(Lvals) # number of Loci values
  num_rho <- length(rhovals) # number of rho values
  # for loop that sets up the different sets of trials (all combinations of Lvals and rhovals)
  for (r in rhovals) for (L in Lvals){
    k <- k+1 
    types[[k]] <- list() 
    types[[k]]$rho = r 
    types[[k]]$numL = L
  }
  # sets up data frame for all z posterior means
  nl <- numeric(0) # will be vector containing values of loci number for all trials
  np <- numeric(0) # will be vector of values for rho values of all trials
  for (L in Lvals){
    a1 <- rep(L,N)
    nl <- c(nl,a1)
  }
  nl <- rep(nl,num_rho)
  for (r in rhovals){
    a2 <- rep(r,num_Ls*N)
    np <- c(np,a2)
  }
  z_df <- data.frame(contam_prob = rep(np,n), loci_number = rep(nl,n), z = rep(0,n*N*num_rho*num_Ls), z_id = rep(0,n*N*num_rho*num_Ls))
  
  # Sets up data frame for the rho values
  rho_v = numeric(0) # vector with all the rho values for the data frame
  for (r in rhovals){
    a3 <- rep(r,num_Ls)
    rho_v <- c(rho_v,a3)
  }
  rho_df <- data.frame(contam_prob = rep(rho_v,n),loci_number = rep(Lvals,n*num_rho), rho_pm = rep(0,n*num_rho*num_Ls))
 
  # Sets up data frame for allele frequencies
  alls <- numeric(0)
  aps <- numeric(0)
  for (L in Lvals){
    a4 <- rep(L,L)
    alls <- c(alls,a4)
  }
  alls <- rep(alls,num_rho)
  s <- sum(Lvals)
  for (r in rhovals){
    a5 <- rep(r,s)
    aps <- c(aps,a5)
  }
  allele_df <- data.frame(rho_value = rep(aps,n), loci_number = rep(alls,n))
  
  total <- N*num_rho*num_Ls
  
  # Runs all MCMCs.  
  # MCMC contains the z posterior means, the rho posterior means, the allele frequencies, the allele estimates, and the upper and lower bound for the 90% allele intervals
  MCMC <- lapply(1:n, function(x) {lapply(types, function(x) test_MCMC(sample_data = sample_data, N = N, l = l, L = x$numL, p = x$rho, inters = inters))})
  alleles <- numeric(0) # correct allele frequencies
  ameans <- numeric(0) # estimated allele frequencies
  top_int <- numeric(0) # upper bound of allele frequencies
  bottom_int <- numeric(0) # lower bound of allele frequencies
  
  # Combines all data so that it can be put into the data frames
  for(i in 1:n){
    values <- MCMC[[i]] # the ith sample MCMC
    #indices for the z to be output
    zindex1 <- (i-1)*total + 1
    zindex2 <- i*total
    
    #indices for the rho values to be output
    pindex1 <- (i-1)*num_rho*num_Ls + 1
    pindex2 <- i*num_rho*num_Ls
    # inputs relevant z values into data frame
    z_df$z[zindex1:zindex2] <- as.vector(sapply(values, function(x) {x$z_pm}))
    z_df$z_id[zindex1:zindex2] <- as.vector(sapply(values, function(x) {x$z_id}))
    # inputs relevant rho values
    rho_df$rho_pm[pindex1:pindex2] <- as.vector(sapply(values, function(x) {x$rho_pm}))
    
    # setting up vectors containing allele data
    a <- unlist(lapply(values, function(x) {x$afreqs}))
    alleles <- c(alleles, a) # true allele frequencies
    am <- unlist(lapply(values, function(x) {x$alle_pm}))
    ameans <- c(ameans, am) # allele frequency estimates
    tp <- unlist(lapply(values, function(x) {x$top_int}))
    top_int <- c(top_int, tp) # upper bound
    bt <- unlist(lapply(values, function(x) {x$bottom_int}))
    bottom_int <- c(bottom_int, bt) # lower bound
  }
  #puts allele data in data frame
  allele_df$alle_freq <- alleles
  allele_df$estimates <- ameans
  allele_df$topint <- top_int
  allele_df$bottomint <- bottom_int
  
  MCMC_zplots(z_df) # plots z values
  MCMC_alleleplot(allele_df) # plots allele data
  MCMC_rhoplot(rho_df) # plots rho data
  allele_table(allele_df,Lvals,rhovals)
  
  list(z = z_df, allele = allele_df, rho = rho_df)
}

  # look into the package parallel => mclapply()
  