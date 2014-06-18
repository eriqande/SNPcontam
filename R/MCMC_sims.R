#' @export
MCMC_sims <- function(sample_data, N, Lvals, rhovals, l, alpha, beta, lambda, inters,n){
  # create list with all simulations
  types <- list() 
  k <- 0
  for (r in rhovals) for (L in Lvals){
    k <- k+1; types[[k]] <- list(); types[[k]]$rho = r; types[[k]]$numL = L
  }
  nl <- rep(c(rep(Lvals[1],N), rep(Lvals[2],N),rep(Lvals[3],N), rep(Lvals[4],N)),4)
  np <- c(rep(rhovals[1],4*N),rep(rhovals[2],4*N),rep(rhovals[3],4*N),rep(rhovals[4],4*N))
  z_df <- data.frame(contam_prob = rep(np,n), loci_number = rep(nl,n), z = rep(0,n*N*16), z_id = rep(0,n*N*16))
  rho_v <- c(rep(rhovals[1],4), rep(rhovals[2],4), rep(rhovals[3],4), rep(rhovals[4],4))
  rho_df <- data.frame(contam_prob = rep(rho_v,n),loci_number = rep(Lvals,4*n), rho_pm = rep(0,n*16))
  alls <- rep(c(rep(Lvals[1],Lvals[1]), rep(Lvals[2],Lvals[2]), rep(Lvals[3],Lvals[3]), rep(Lvals[4],Lvals[4])),4)
  aps <- c(rep(rhovals[1],sum(Lvals)), rep(rhovals[2],sum(Lvals)), rep(rhovals[3],sum(Lvals)), rep(rhovals[4],sum(Lvals)))
  allele_df <- data.frame(rho_value = rep(aps,n), loci_number = rep(alls,n), alle_freq = rep(0,4*sum(Lvals)*n), estimates = rep(0,4*sum(Lvals)*n))
  total <- N*16
  MCMC <- lapply(1:n, function(x) {lapply(types, function(x) test_MCMC(sample_data = sample_data, N = N, l = l, L = x$numL, p = x$rho, inters = 1000))})
  alleles <- numeric(0)
  ameans <- numeric(0)
  top_int <- numeric(0)
  bottom_int <- numeric(0)
  for(i in 1:n){
    values <- MCMC[[i]]
    zindex1 <- (i-1)*total + 1
    zindex2 <- i*total
    pindex1 <- (i-1)*16 + 1
    pindex2 <- i*16
    z_df$z[zindex1:zindex2] <- as.vector(sapply(values, function(x) {x$z_pm}))
    z_df$z_id[zindex1:zindex2] <- as.vector(sapply(values, function(x) {x$z_id}))
    rho_df$rho_pm[pindex1:pindex2] <- as.vector(sapply(values, function(x) {x$rho_pm}))
    a <- unlist(lapply(values, function(x) {x$afreqs}))
    alleles <- c(alleles, a)
    am <- unlist(lapply(values, function(x) {x$alle_pm}))
    ameans <- c(ameans, am)
    tp <- unlist(lapply(values, function(x) {x$top_int}))
    top_int <- c(top_int, tp)
    bt <- unlist(lapply(values, function(x) {x$bottom_int}))
    bottom_int <- c(bottom_int, bt)
  }
  allele_df$alle_freq <- alleles
  allele_df$estimates <- ameans
  allele_df$topint <- top_int
  allele_df$bottomint <- bottom_int
  
  MCMC_zplots(z_df)
  MCMC_alleleplot(allele_df)
  MCMC_rhoplot(rho_df)
  
  list(z = z_df, allele = allele_df, rho = rho_df)
}

  # look into the package parallel => mclapply()
  