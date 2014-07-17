

#' read in standard two-column SNP data and spit back allele freqencies
#' 
#' @export
allele_frequencies_from_two_column_data <- function(S) {
  # make a list of unique alleles
  uniq_alleles <- lapply(seq(1, ncol(S), 2), function(x) levels(factor(c(S[,x], S[, x+1]))))
  names(uniq_alleles) <- colnames(S)[seq(1, ncol(S), 2)]
  
  # drop loci that do not have two alleles
  Have2 <- sapply(uniq_alleles, function(x) length(x)==2)
  Have2rep <- rep(Have2, each=2)  # this can be used to extract the loci from M2 and B2 that have two alleles
  
  DropTheseLoci <- names(uniq_alleles)[!Have2]
  if(length(DropTheseLoci)>0) warning(paste("Dropping these loci that have > of < then 2 alleles:", paste(DropTheseLoci, collapse=", ")))
  
  S2 <- S[, Have2rep]
  uniq_alleles2 <- uniq_alleles[Have2]
  
  # now, make sure that the levels of S2's alleles are appropriately set
  for(i in seq(1, ncol(S2), 2)) {
    S2[, i] <- factor(S2[, i], levels = uniq_alleles2[[colnames(S2)[i]]])
    S2[, i+1] <- factor(S2[, i+1], levels = uniq_alleles2[[colnames(S2)[i]]])
  }
  
  
  # now count the number of zero and one alleles in each population at each locus in baseline B3
  alle_counts <- lapply(seq(1, ncol(S2), 2), function(y) table(c(S2[, y], S2[, y+1])))
  names(alle_counts) <- colnames(S2)[seq(1, ncol(S2), 2)]
  alle_counts_array <- simplify2array(alle_counts)
  
  
  totsums <- colSums(alle_counts_array)
  
  proportions <- t(apply(alle_counts_array, 1, function(x) x/totsums))
  dimnames(proportions) <- list(Alleles = c("0", "1"), Loci = colnames(proportions))
  
}