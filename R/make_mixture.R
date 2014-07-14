#add variable for index of fish to put in mixture (clean and contaminated)
# use colnames of noncontam
make_mixture <- function(baseline, N, p, clean_list = NULL, contam_list = NULL){
  N_c = N*p
  swfs <- baseline
  rownames(swfs) <- swfs$ID
  swfs2 <- swfs
  swfs2$Pop <- factor(as.character(swfs2$Pop), levels = unique(as.character(swfs2$Pop)))
  swfs2$RepPop <- factor(as.character(swfs2$RepPop), levels = unique(as.character(swfs2$RepPop)))
  
#### Pull out Contaminated Genotypes and make a mixture sample with them ####
  L <- ncol(swfs2)-4
  total <- dim(swfs2)[1]
  cc <- 1:nrow(swfs2); names(cc) <- swfs2$ID 
  if (p != 0){
    #gets contaminted sample indices randomly if there is no input in contam_list
    if (length(contam_list) == 0){
      contam <- sample(1:total,2 * N_c)
      contam1 <- contam[1:N_c]
      contam2 <- contam[(N_c+1):(2*N_c)]
    # uses the list of indices to get contaminated samples
    }else{
      contam1 <- unname(cc[contam_list[1,]])
      contam2 <- unname(cc[contam_list[2,]])
      contam <- c(contam1,contam2)
      }
    # gets clean sample indices randomly if there is no clean_list
    if(length(clean_list) == 0){
      a <- sample((1:total)[-contam], N - N_c)
    # gets the list of indices for clean samples
    }else{
      a <- unname(cc[clean_list])
    }
    
    genos <- swfs2[a, -4]  # uncontaminated fish in the mixture
    geno1 <- swfs2[contam1,-(1:4)]  # for contaminated fish in mixture
    geno2 <- swfs2[contam2,-(1:4)]
    g1_info <- swfs2[contam1,(1:3)]
    g2_info <- swfs2[contam2,(1:3)]
    
    contaminate_genos <- function(x,y){
      geno_c <- x # geno_c will be the final version of the contaminated genotypes
      alle1 <- seq(1,(L-1),2) # odd rows which correspond with first gene copy
      alle2 <- seq(2,L,2) # even rows which correspond with second gene copy
      dif1 <- which((x[,alle1] != y[,alle1]) == TRUE, arr.ind = TRUE) # indices of first gene copies that are different between individuals
      dif2 <- which((x[,alle2] != y[,alle2] & x[,alle1] == y[,alle1]) == TRUE, arr.ind = TRUE) # indices of 2nd gene copies different between individuals 
      missing <- which(is.na(y) == TRUE, arr.ind=TRUE) # any missing values from the second list of genotypes
      #missing <- which(is.na(x) == TRUE, arr.ind=TRUE)
      # Code below creates the contamination
      for(i in 1:nrow(dif1)){geno_c[(dif1[i,1]),(2*dif1[i,2])] = y[dif1[i,1],(2*dif1[i,2]-1)]}
      for(j in 1:nrow(dif2)){geno_c[(dif2[j,1]),(2*dif2[j,2]-1)] = y[dif2[j,1],(2*dif2[j,2])]}
      if(length(missing) != 0){
        for(k in 1:nrow(missing)){geno_c[(missing[k,1]),(missing[k,2])] = NA}
        # for(k in 1:nrow(missing)){geno_c[missing[k,1],missing[k,2]] = y[missing[k,1],missing[k,2]]}
      }
      geno_c
      }
    
    geno_c0 <- contaminate_genos(geno1,geno2) # contaminated samples
    rownames(geno_c0) <- paste(rownames(geno1), rownames(geno2), sep="-x-")
    info <- sapply(1:3, function(x) paste(g1_info[,x], g2_info[,x], sep = "-x-"))
    geno_c1 <- cbind(info,geno_c0)
    colnames(geno_c1)[1:3] = c("RepPop","RepUnit","Pop")
    
    data <- rbind(genos,geno_c1) # mixture
  
  # in the end, we drop those from the baseline too
    bline <- swfs2[-c(a, contam),-4]
    }

#### Make Mixture and Baseline if There is No Contamination ####
  else {
    if (length(clean_list) == 0){
      a <- sample((1:total), N)
    }else{
      cc <- 1:nrow(swfs2); names(cc) <- swfs2$ID; a <- unname(cc[clean_list])
    }
    genos <- swfs2[a,-(1:4)]
    snp_genos <- get_snp_genos(genos)
    data <- snp_genos$mat
    bline <- swfs2[-c(a),-4]
    
  }

#### Output ####
list(bline = bline, mixture = data)
}