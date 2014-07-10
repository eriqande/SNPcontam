#add variable for index of fish to put in mixture (clean and contaminated)
# use colnames of noncontam
make_mixture <- function(baseline, N, p){
  N_c = N*p
  swfs <- baseline
  rownames(swfs) <- swfs$ID
  swfs2 <- swfs
  swfs2$Pop <- factor(as.character(swfs2$Pop), levels = unique(as.character(swfs2$Pop)))
  swfs2$RepPop <- factor(as.character(swfs2$RepPop), levels = unique(as.character(swfs2$RepPop)))
  
#### Pull out Contaminated Genotypes and make a mixture sample with them ####
  total <- dim(swfs2)[1]
  if (p != 0){
  contam <- sample(1:total,2 * N_c)
  contam1 <- contam[1:N_c]
  contam2 <- contam[(N_c+1):(2*N_c)]
  
  a <- sample((1:total)[-contam], N - N_c)
  
  genos <- swfs2[a,-(1:4)]   # uncontaminated fish in the mixture
  geno1 <- get_snp_genos(swfs2[contam1,-(1:4)])$mat  # for contaminated fish in mixture
  geno2 <- get_snp_genos(swfs2[contam2,-(1:4)])$mat
  geno_c <- geno1 + geno2
  geno_c[geno_c == 4] = 2
  geno_c[geno_c == 3] = 1
  geno_c[geno1 ==1 & geno2 == 1] = 1
  geno_c[(geno1 == 0 & geno2 == 2)|(geno1 == 2 & geno2 == 0)] = 1
  geno_c[is.na(geno1) & is.na(geno2) == FALSE] = geno2[is.na(geno1) & is.na(geno2) == FALSE]
  geno_c[is.na(geno2) & is.na(geno1) == FALSE] = geno1[is.na(geno2) & is.na(geno1) == FALSE]
  
  colnames(geno_c) <- paste(colnames(geno1), colnames(geno2), sep="-x-")
  
  snp_genos <- get_snp_genos(genos)
  data <- cbind(snp_genos$mat,geno_c)
  
  # in the end, we drop those from the baseline too
  bline <- swfs2[-c(a, contam),-4]
  }

#### Make Mixture and Baseline if There is No Contamination ####
  else {
  a <- sample((1:total), N)
  genos <- swfs2[a,-(1:4)]
  snp_genos <- get_snp_genos(genos)
  data <- snp_genos$mat
  bline <- swfs2[-c(a),-4]
  }

#### Output ####
list(bline = bline, mixture = data)
}