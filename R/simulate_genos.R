#' @export
simulate_genos <- function(N,L,p,l, afreqs){
  n_contam <- round(N*p,0)
  l_contam <- round(L*l,0)
  g <- matrix(0,L,N)
  # get genotype frequencies
  all_L <- ncol(afreqs)
  loci <- sample(1:all_L,L,replace = TRUE)
  genes <- afreqs[2,loci] + runif(L,-.01,.01)
  genes[genes==0.0] <- .01
  genes[genes<0.0] <- -genes[genes<0.0]
  genes[genes>1.0] <- 2.0 - genes[genes>1.0]
  gfreqs <- likelihood(genes)

  # randomly pick individuals to be contaminated
  if (p==0){
    for (i in 1:L){
      g[i,] <- sample(0:2,N,TRUE,prob=gfreqs$clean_prob[,i])
    }
  } else{
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
  }
  if (p==0){
  contam_id <- 0
  } else{
    contam_id = id
  }
  list(geno = g, contam_id = contam_id, loci = loci, afreqs = genes)
}