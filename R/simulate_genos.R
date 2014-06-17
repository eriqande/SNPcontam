#' @export
simulate_genos <- function(N,L,p,l,sample_data){
  n_contam <- round(N*p,0)
  l_contam <- round(L*l,0)
  g <- matrix(0,L,N)
  # get data and allele frequencies
  snp_genos <- get_snp_genos(sample_data)
  snp_indices <- genos_to_indicators(g = snp_genos$mat)
  geno_counts <- count_genos(snp_indices)
  afreqs <- alle_freqs(geno_counts)
  # get genotype frequencies
  all_L <- nrow(snp_genos$mat)
  loci <- sample(1:all_L,L,replace = TRUE)
  genes <- afreqs[2,loci]
  gfreqs <- likelihood(genes)
  
  # randomly pick individuals to be contaminated
  id <- sample(1:N,n_contam,replace = FALSE)
  contam_loci <- sample(1:L,l_contam,replace = FALSE)
  
  #simulate genotypes
  for (i in 1:L){
    if (i %in% contam_loci){
    # stores genotype at locus in all individuals
    g[i,id] <- sample(0:2, n_contam, TRUE, prob = gfreqs$contam_prob[,i])
    g[i,-(id)] <- sample(0:2,(N-n_contam),TRUE,prob = gfreqs$clean_prob[,i])
  } else{
    g[i,] <- sample(0:2,N,TRUE,prob=gfreqs$clean_prob[,i])
  }
  }
  list(geno = g, contam_id = id, loci = loci, afreqs = genes)
}