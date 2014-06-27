MCMC_sims <- function(sample_data, N, Lvals, rhovals, l=1, alpha=0.5, beta=0.5, lambda=0.5, inters=1000,n){
  snp_genos <- get_snp_genos(sample_data)
  snp_indices <- genos_to_indicators(g = snp_genos$mat)
  geno_counts <- count_genos(snp_indices)
  afreqs <- alle_freqs(geno_counts)
  
  # create list with all simulations
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
  
  # Runs all MCMCs.  
  # MCMC contains the z posterior means, the rho posterior means, the allele frequencies, the allele estimates, and the upper and lower bound for the 90% allele intervals
  MCMC <- mclapply(types, function(x) {
      lapply(1:n, function(y) {
          list(params = x, output = test_MCMC(afreqs = afreqs, N = N, l = l, L = x$numL, p = x$rho, inters = inters))
          })
      },mc.cores=length(types))
 
  # Makes the z output data frame
  slurp_mcmc_z_output <- function(y) {
    data.frame(z = y$output$z_pm, z_id = y$output$z_id, contam_prob = y$params$rho, loci_number = y$params$numL )
  }
  ztmp1 <- mclapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x) {df <- slurp_mcmc_z_output(x); df$rep_num=rep; df }) }) 
  ztmp2 <- unlist(ztmp1, recursive = FALSE)
  z_df <- do.call(what = rbind, args = ztmp2)
  
  # Makes the rho output data frame
  slurp_mcmc_rho_output <- function(y){
    data.frame(rho_pm = y$output$rho_pm, contam_prob = y$params$rho, loci_number = y$params$numL)
  }
  rtmp1 <- mclapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x) {df <- slurp_mcmc_rho_output(x); df$rep_num=rep; df})})
  rtmp2 <- unlist(rtmp1, recursive = FALSE)
  rho_df <- do.call(what = rbind, args = rtmp2)
  
  # Makes the allele output data frame
  slurp_mcmc_allele_output <- function(y){
    data.frame(alle_freq = y$output$afreqs, estimates = y$output$alle_pm, topint = y$output$top_int, bottomint = y$output$bottom_int, contam_prob = y$params$rho, loci_number = y$params$numL)
  }
  atmp1 <- mclapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x){df <- slurp_mcmc_allele_output(x); df$rep_num=rep; df})})
  atmp2 <- unlist(atmp1, recursive = FALSE)
  allele_df <- do.call(what = rbind, args = atmp2)
    
  list(types = types, z = z_df, allele = allele_df, rho = rho_df, rhovals = rhovals, Lvals = Lvals)
} 