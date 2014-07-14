# Make the Baseline and the Mixture
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

# Get clean_list and contam list for simulation
mixed_lists <- function(baseline, props, N ,p , less.NA = NULL){
  if(length(less.NA) > 0){
      num_na <- rowSums(is.na(baseline))
      newb <- baseline[-which(num_na > less.NA),]
  }
  tmp <- sapply(1:length(props), function(x) {t_pop <- sum(newb$RepPop == levels(newb$RepPop)[x]); f_prop <- props[x]/t_pop; pop_props <- rep(f_prop,t_pop); pop_props})
  bline_probs <- matrix(unlist(tmp),ncol = 1)
  rownames(bline_probs) <- newb$ID
  contam_list <- matrix(sample(x = newb$ID, size = 2*N*p, replace = TRUE, prob = bline_probs),nrow = 2)
  bline_probs2 <- bline_probs[-which(rownames(bline_probs) %in% contam_list)]
  clean_list <- sample(x = newb$ID[-which(newb$ID %in% contam_list)], size = N-N_c, replace = TRUE,prob = bline_probs2)
  list(clean_list = clean_list, contam_list = contam_list)
}

# Simulation Function
mixed_MCMC_sims <- function(baseline, N, p, fish_pops,inters,contamination = TRUE, less.NA = NULL, n){
  
  # function for simplifying MCMC data
  mixed_MCMC_with_means <- function(baseline, N, p, contamination = TRUE, inters, props, less.NA = NULL, alpha=.5, beta=.5){
    N_c = N*p
    mixture_lists <- mixed_lists(baseline, props, N, p, less.NA = less.NA)
    b <- make_mixture(baseline,N,p,clean_list = mixture_lists$clean_list, contam_list = mixture_lists$contam_list)
    MCMC_mats <- prepare_base_and_mix_for_mixed_MCMC(B = b$bline, B_locstart = 4, B_pops = b$bline$Pop, M = b$mixture, M_locstart = 4)
    contam_mat <- Pcontam(snp_zeroes = MCMC_mats$zero, snp_ones = MCMC_mats$ones, genos = MCMC_mats$mixmat, lambda = .5)
    clean_mat <- P_likelihood(snp_zeroes = MCMC_mats$zero, snp_ones = MCMC_mats$ones, genos = MCMC_mats$mixmat, lambda = .5)
    MCMC_data <- mixed_MCMC(MCMC_mats$mixmat, contam_mat, clean_mat, inters = inters, contamination = TRUE)
    z_pm <- colMeans(MCMC_data$z)
    z_id <- c(rep(0,N),rep(1,N*p))
    rho_pm <- mean(MCMC_data$prob_contam)
    P = nrow(clean_mat)
    # population assignment for clean populations
    if (contamination){
      clean_pops <- lapply(1:(N-N_c), function(x) {MCMC_data$pops[,x][2*which(MCMC_data$z[,x] == 0) - 1]})
      tmp_pops <- lapply(clean_pops, function(x) {factor(x,level = 1:P)}) 
      freq <- lapply(tmp_pops, function(x) {as.numeric(table(x))})
      pop_means <- do.call(what = rbind, args = lapply(freq, function(x) x))/inters
      colnames(pop_means) <- levels(b$bline$RepPop)
      if (p == 0) {rownames(pop_means) <- rownames(b$mixture)}else {rownames(pop_means) <- rownames(b$mixture)[1:(N-N_c)]}
    }else {
      tmp_pops <- lapply(1:N, function(x) {factor(MCMC_data$pop[,x],level = 1:P)}) 
      freq <- lapply(tmp_pops, function(x) {as.numeric(table(x))})
      pop_means <- do.call(what = rbind, args = lapply(freq, function(x) x))/inters
      colnames(pop_means) <- levels(b$bline$RepPop)
      rownames(pop_means) <- rownames(b$mixture)
    }
    
    #contaminated population data
    contam_pops <- lapply((N-N_c + 1):N, function(x) {if (z_pm[x] == 1){MCMC_data$pop[,x]} else{MCMC_data$pops[,x][-c((2*which(MCMC_data$z[,x] == 0) - 1),(2*which(MCMC_data$z[,x] == 0)))]}})
    
    list(z_pm = z_pm , z_id = z_id, rho_pm = rho_pm, clean_means = pop_means, contam_pops = contam_pops, MCMC = MCMC_data, mixture_ids = rownames(b$mixture))  
  }
  
  # run MCMC on all simulation mixtures
  types <- lapply(1:nrow(fish_pops), function(x) list(fish_pops[x,]))
  names(types) <- rownames(fish_pops)
  MCMC <- lapply(types, function(x) {
    lapply(1:n, function(y) {
      list(params = names(x), output = mixed_MCMC_with_means(baseline = baseline, N = N, p = p, contamination = contamination, inters = inters, props = unlist(x), less.NA = less.NA))
    })
  }) #mc.cores=length(types))
  
  # MAkes the z output data frame
  slurp_mcmc_z_output2 <- function(y) {
    data.frame(z = y$output$z_pm, z_id = y$output$z_id, fishery_props = y$params)
  }
  ztmp1 <- mclapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x) {df <- slurp_mcmc_z_output(x); df$rep_num=rep; df }) }) 
  ztmp2 <- unlist(ztmp1, recursive = FALSE)
  z_df <- do.call(what = rbind, args = ztmp2)
  
  # Makes the rho output data frame
  slurp_mcmc_rho_output2 <- function(y){
    data.frame(rho_pm = y$output$rho_pm, cfishery_props = y$params)
  }
  rtmp1 <- mclapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x) {df <- slurp_mcmc_rho_output(x); df$rep_num=rep; df})})
  rtmp2 <- unlist(rtmp1, recursive = FALSE)
  rho_df <- do.call(what = rbind, args = rtmp2)
  
  list(z_df = z_df, rho_df = rho_df)
}