get_likelihood_matrices <- function(bline, mixture){
  get_zeroes_and_ones <- function(x) {
    y <- x[, -(1:2)]
    snp_genos <- get_snp_genos(y)$mat
    snp_indics <- genos_to_indicators(g = snp_genos)
    geno_counts <- count_genos(snp_indics)
    af <- t(alle_freqs(geno_counts, proportion = F))
  }
  
  # Get zero and one allele counts for each population
  pop.list <- split(bline, bline$Pop)
  alle.counts.list <- lapply(pop.list, get_zeroes_and_ones)
  zeros <- do.call(what = cbind, args = lapply(alle.counts.list, function(x) x[,"0"]))
  ones <-  do.call(what = cbind, args = lapply(alle.counts.list, function(x) x[,"1"]))
  
  
  # Getting the two probability matrices ####
  clean_prob <- P_likelihood(zeros,ones,mixture,.5)
  contam_prob <- Pcontam(zeros,ones,mixture,.5)
  
  list(clean_prob = clean_prob, contam_prob = contam_prob)
}