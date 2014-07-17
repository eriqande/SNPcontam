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
    data <- swfs2[a,-4]
    bline <- swfs2[-c(a),-4]
    
  }
  
  #### Output ####
  list(bline = bline, mixture = data)
}

#### Get clean_list and contam list for simulation ####
mixed_lists <- function(baseline, props, N ,p , less.NA = NULL){
  N_c = N*p
  if(length(less.NA) > 0){
      num_na <- rowSums(is.na(baseline))
      newb <- baseline[-which(num_na > less.NA),]
  }
  tmp <- sapply(1:length(props), function(x) {t_pop <- sum(newb$RepPop == levels(newb$RepPop)[x]); f_prop <- props[x]/t_pop; pop_props <- rep(f_prop,t_pop); pop_props})
  bline_probs <- matrix(unlist(tmp),ncol = 1)
  rownames(bline_probs) <- newb$ID
  if(p!=0) {
    contam_list <- matrix(sample(x = newb$ID, size = 2*N*p, replace = TRUE, prob = bline_probs),nrow = 2)
    bline_probs2 <- bline_probs[-which(rownames(bline_probs) %in% contam_list)]
    clean_list <- sample(x = newb$ID[-which(newb$ID %in% contam_list)], size = N-N_c, replace = TRUE,prob = bline_probs2)
  } else{
    clean_list <- sample(x = newb$ID, size = N, replace = TRUE, prob = bline_probs)
    contam_list = NULL
  }
  list(clean_list = clean_list, contam_list = contam_list)
}

#### Simulation Function ####
mixed_MCMC_sims <- function(baseline, N, p, fish_pops,inters,contamination = TRUE, less.NA = NULL, n, MAX_CORES = 4){
  P = length(levels(baseline$Pop))
  # function for simplifying MCMC data
  mixed_MCMC_with_means <- function(baseline, N, p, contamination = TRUE, inters, props, less.NA = NULL, alpha=.5, beta=.5, burnin = 100){
    N_c = N*p
    mixture_lists <- mixed_lists(baseline, props, N, p, less.NA = less.NA) # get lists of samples to using in mixture
    
    # get baseline and mixture using the sampled lists from mixture_lists
    b <- make_mixture(baseline,N,p,clean_list = mixture_lists$clean_list, contam_list = mixture_lists$contam_list)
    
    # prepare matrices for the MCMC
    MCMC_mats <- prepare_base_and_mix_for_mixed_MCMC(B = b$bline, B_locstart = 4, B_pops = b$bline$Pop, M = b$mixture, M_locstart = 4)
    
    # get clean and contam likelihood matrices and run MCMC
    contam_mat <- Pcontam(snp_zeroes = MCMC_mats$zero, snp_ones = MCMC_mats$ones, genos = MCMC_mats$mixmat, lambda = .5)
    clean_mat <- P_likelihood(snp_zeroes = MCMC_mats$zero, snp_ones = MCMC_mats$ones, genos = MCMC_mats$mixmat, lambda = .5)
    MCMC_data <- mixed_MCMC(MCMC_mats$mixmat, contam_mat, clean_mat, inters = inters, contamination = TRUE)
    
    # get posterior means for z and rho
    z_pm <- colMeans(MCMC_data$z[-(1:burnin),])
    z_id <- c(rep(0,(N-N_c)),rep(1,N*p))
    rho_pm <- mean(MCMC_data$prob_contam[-(1:burnin)])
    P = nrow(clean_mat)
    
    # Posterior means for mixing proportions (pii)
    pii = colMeans(MCMC_data$mixing[-(1:burnin),])
    
    # Population assignment for clean populations
    clean_pops <- lapply(1:(N-N_c), function(x) {MCMC_data$pops[-(1:(2*burnin)),x][2*which(MCMC_data$z[-(1:(2*burnin)),x] == 0) - 1]})
    tmp_pops <- lapply(clean_pops, function(x) {factor(x,level = 1:P)}) 
    freq <- lapply(tmp_pops, function(x) {as.numeric(table(x))})
    # fraction of assignments of each individual to each population
    pop_means <- do.call(what = rbind, args = lapply(freq, function(x) x))/inters
    colnames(pop_means) <- levels(b$bline$RepPop)
    if (p == 0) {rownames(pop_means) <- rownames(b$mixture)}else {rownames(pop_means) <- rownames(b$mixture)[1:(N-N_c)]}
    
    # contaminated population data
    if (p != 0) {
      contam_pops <- lapply((N-N_c + 1):N, function(x) {if (z_pm[x] == 1) {MCMC_data$pop[-(1:(2*burnin)),x]} else{MCMC_data$pops[-(1:(2*burnin)),x][-c((2*which(MCMC_data$z[-(1:(2*burnin)),x] == 0) - 1),(2*which(MCMC_data$z[-(1:(2*burnin)),x] == 0)))]}})
    
    # Contam population assignment
      opts <- combinations(P,2,1:P, repeats.allowed = TRUE) # all combinations of 2 populations
      t_opts <- nrow(opts) # number of combinations
      # Matrix of each individuals contamination population data. Two columns that are the two different combinations.
      mat_ind <- lapply(contam_pops, function(x) matrix(x,ncol=2, byrow = TRUE))
      # Number of times each of the combinations appears
      counts <- lapply(mat_ind, function(y) {sapply(1:t_opts, function(x) {sum(y[,1] == opts[x,1] & y[,2] == opts[x,2] | y[,1] == opts[x,2] & y[,2] == opts[x,1])})})
      # Pair with the maximum assignments for each individual
      max_pair <- sapply(counts, function(x) {a <- which(max(x) == x); if(length(a) > 1){a2 <- a[[1]]} 
                                                else{a2 <- a}; opts[a2,]})
      # Counts without the maximum
      counts2 <- lapply(counts, function(x) {a <- which(max(x) == x); if(length(a) > 1){a2 <- a[[1]]} 
                                               else{a2 <- a}; x[-a2]})
      # Pair with the second highest counts
      max_pair2 <- sapply(counts2, function(x) {a <- which(max(x) == x); if(length(a) > 1){a2 <- a[[1]]} 
                                                  else{a2 <- a}; opts[a2,]})
      # Population Natmes
      populations = levels(b$bline$RepPop)
      #Data frame wiht all of the name of the individual and the populations of origin according to the MCMC
      contam_df <- data.frame(max_pair = matrix(populations[t(max_pair)],ncol=2), 
                   max_pair2 = matrix(populations[t(max_pair2)], ncol=2))
    } else{contam_pops = NULL; contam_df = NULL}
    
    list(z_pm = z_pm , z_id = z_id, rho_pm = rho_pm, clean_means = pop_means, contam_pops = contam_pops, contam_df = contam_df, MCMC = MCMC_data, mixing_props = pii, mixture_ids = rownames(b$mixture))  
  }
  
  # run MCMC on all simulation mixtures
  # types is all of the proportion types 
  f_types <- lapply(1:nrow(fish_pops), function(x) list(pop = fish_pops[x,], name = rownames(fish_pops)[x]))
  types <- list() 
  k <- 0
  for (r in p) for (f in f_types){
    k <- k+1 
    types[[k]] <- list() 
    types[[k]]$rho = r   
    types[[k]]$f_types = f
  }
  MCMC <- mclapply(types, function(x) {
    lapply(1:n, function(y) {
      list(params = list(name = x$f_types$name, proportions = x$f_types$pop, rho = x$rho), 
           output = mixed_MCMC_with_means(baseline = baseline, N = N, p = x$rho, contamination = contamination, inters = inters, props = unlist(x$f_types$pop), less.NA = less.NA))
    })
  }, mc.cores=min(length(types), MAX_CORES))
  
  save(MCMC, file="MCMC.rda")
  
  # MAkes the z output data frame
  slurp_mcmc_z_output2 <- function(y) {
    data.frame(z = y$output$z_pm, z_id = y$output$z_id, fishery = y$params$name, ID = y$output$mixture_ids, rho = y$params$rho, stringsAsFactors = FALSE)
  }
  ztmp1 <- mclapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x) {df <- slurp_mcmc_z_output2(x); df$rep_num=rep; df })}, mc.cores = min(length(MCMC), MAX_CORES)) 
  ztmp2 <- unlist(ztmp1, recursive = FALSE)
  z_df <- do.call(what = rbind, args = as.list(ztmp2))
  
  # Makes the rho output data frame
  slurp_mcmc_rho_output2 <- function(y){
    data.frame(rho_pm = y$output$rho_pm, rho = y$params$rho, fishery = y$params$name, stringsAsFactors = FALSE)
  }
  rtmp1 <- mclapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x) {df <- slurp_mcmc_rho_output2(x); df$rep_num=rep; df})}, mc.cores = min(length(MCMC), MAX_CORES))
  rtmp2 <- unlist(rtmp1, recursive = FALSE)
  rho_df <- do.call(what = rbind, args = as.list(rtmp2))
  
  # MAkes the mixing output data frame
  slurp_mcmc_pii_output <- function(y){
    colnames(y$params$proportions) <- NULL
    rownames(y$params$proportions) <- NULL
    data.frame(populations = colnames(y$output$clean_means),true_mixing = t(y$params$proportions), mixing_pm = y$output$mixing_props, rho = y$params$rho, fishery = y$params$name, stringsAsFactors = FALSE)
    }
  mtmp1 <- mclapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x) {df <- slurp_mcmc_pii_output(x); df$rep_num=rep; df})}, mc.cores = min(length(MCMC), MAX_CORES))
  mtmp2 <- unlist(mtmp1, recursive = FALSE)
  pii_df <- do.call(what = rbind, args = as.list(mtmp2))
    
  # Makes data frame for clean population means
  slurp_mcmc_u_output <- function(y){
    data.frame(fishery = y$params$name, u = y$output$clean_means, rho = y$params$rho, stringsAsFactors = FALSE)
    }
#  utmp1 <- mclapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x) {df <- slurp_mcmc_u_output(x); df$rep_num=rep; df})}, mc.cores = min(length(MCMC), MAX_CORES))
  utmp1 <- lapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x) {df <- slurp_mcmc_u_output(x); df$rep_num=rep; df})})
  utmp2 <- unlist(utmp1, recursive = FALSE)
  u_df <- do.call(what = rbind, args = as.list(utmp2))
    
  #Contam Data
  #Find maximum pair
  slurp_mcmc_u_contam_output <- function(y,N){
    if(y$params$rho != 0){
      p = y$params$rho
      N_c = N*p
      data.frame(ID = y$output$mixture_ids[(N-N_c+1):N], fishery = y$params$name, u = y$output$contam_df, rho = y$params$rho, stringsAsFactors = FALSE)
    } else{NULL}
  }
  #pair_tmp <- mclapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x) {slurp_mcmc_u_contam_output(x,N)})}, mc.cores = min(length(MCMC), MAX_CORES))
  pair_tmp <- lapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x) {slurp_mcmc_u_contam_output(x,N)})})
  pair_tmp2 <- unlist(pair_tmp, recursive = FALSE)
  contam_u_df <- do.call(what = rbind, args = as.list(pair_tmp2))
  
  #Just get all of the contam data
#  contam_data <- mclapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x) {x$output$contam_pops})}, mc.cores = min(length(MCMC), MAX_CORES))  
contam_data <- lapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x) {x$output$contam_pops})})  
  
  list(rhovals = p, z_df = z_df, rho_df = rho_df, pii_df, u_df = u_df , contam_u_df = contam_u_df, contam_data = contam_data)
}

mixed_MCMC_rhoplot <- function(rho_df,rhovals,outpath, width = 5, height = 3){
  a <- ggplot(rho_df, aes(x=factor(fishery), y=rho_pm)) + geom_boxplot(outlier.shape=1, outlier.size=1.5) +
    facet_grid(.~rho) + ylab("Posterior Mean") + xlab("Fishery") + geom_hline(aes(yintercept=rho),lty="dotdash") +
    theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.title.x=element_text(size=10), plot.title = element_text(hjust = 0))
  pdf(outpath,width=5,height=3)
  print(a)
  dev.off()
}

mixed_MCMC_population_info <- function(u_df){
  get_pop_data <- function(split_df){
    nc <- ncol(split_df)
    rho <- split_df$rho[1]
    fishery <- as.character(split_df$fishery[1])
    N <- nrow(split_df)
    pop_means <- split_df[,-c(1,(nc-1),nc)]
    tmp_id <- sapply(1:N, function(x) colnames(pop_means)[which(max(pop_means[x,]) == pop_means[x,])])
      
    mc_unit <- strsplit(as.character(tmp_id),("\\."))
    mc_runit <- lapply(mc_unit, function(x) x[2])
    mc_pop <- lapply(mc_unit, function(x) x[4])
    mix <- strsplit(rownames(split_df), "[--:]")
    mix_runit <- lapply(mix, function(x) x[1])
    mix_pop <- lapply(mix, function(x) x[3])
    correct_Pop = sum(unlist(mix_pop) == unlist(mc_pop))/N
    correct_Rep = sum(unlist(mix_runit) == unlist(mc_runit))/N
    list(correct_Pop = correct_Pop, correct_Rep = correct_Rep, rho = rho, fishery = fishery)
  }
  factor(u_df$rho)
  factor(u_df$fishery)
  tmp_pops <- split(u_df, list(u_df$rho, u_df$fishery))
  pop_ids <- lapply(tmp_pops, function(x) get_pop_data(split_df = x))
  
  tmp_pops2 <- lapply(pop_ids, function(x) unlist(x))
  pop_data <- do.call(what = rbind, args = as.list(tmp_pops2))
  pop_df <- data.frame(fishery = pop_data[,4], rho = as.numeric(pop_data[,3]), correct_Pop = as.numeric(pop_data[,1]), correct_RepUnit = as.numeric(pop_data[,2]))
  rownames(pop_df) = NULL
  pop_df
}

mixed_MCMC_ztable <- function(z_df, PPlim = 0.5, outpath) {
  zdata <- function(z, PL = 0.5){
  # NP <== "No Posterior"  (We do not include posterior >0.9 as a factor)
  NP <- table(list(rho=z$rho, L=z$fishery, Z=z$z_id))
  
  # WP <== "with posterior >PPlim" included as a factor
  WP <- table(list(rho=z$rho, L=z$fishery, Z=z$z_id, GreaterThanPPlim=z$z>PL) )
  
  # uncontams with post prob > PPlim
  FalsePos <- WP[,,"0","TRUE"] / NP[,,"0"]
  
  # true contams with post prob > PPlim
  TruePos <- WP[,,"1","TRUE"] / NP[,,"1"]
  
  list(FP = FalsePos, TP = TruePos)}
  tmp_fishery <- strsplit(z_df$fishery,"_")
  z_df$fishery <- sapply(tmp_fishery, function(x) unlist(paste(x[1],x[2],sep = ".")))
  z_tab <- zdata(z = z_df, PL = PPlim)
  ztab_file <- outpath
  cat("{\\bf (a)} Contaminated samples \n", file = ztab_file)
  cat("\\begin{center}\n", file = ztab_file, append = T)
  make_tabular(x = z_tab$TP[-1,], dd = 3, leftcolhead = "$\\rho~~~~~~~$", rightcolshead = "Fishery", outfile = ztab_file, append = T)
  cat("\\end{center}\n", file = ztab_file, append = T)
  cat("{\\bf (b)} Non-contaminated samples \n", file = ztab_file, append = T)
  cat("\\begin{center}\n", file = ztab_file, append = T)
  make_tabular(x = z_tab$FP, dd = 4, leftcolhead = "$\\rho~~~~~~~$", rightcolshead = "Fishery", outfile = ztab_file, append = T)
  cat("\\end{center}\n", file = ztab_file, append = T)
}

make_tabular <- function(x, dd=2, leftcolhead, rightcolshead, outfile="", open_tabular = TRUE, close_tabular = TRUE, append = FALSE) {
  nr <- nrow(x)
  nc <- ncol(x)
  xf <- format(round(x, digits=dd), digits=dd)  # make them strings formatted as desired
  xf <- cbind(rownames(xf), xf)  # add the rownames as a column
  xf <- rbind(c(leftcolhead, colnames(x)), xf) # add the colnames as the top row
  
  
  # remove outfile if it exists
  if(append == FALSE) {if(file.exists(outfile)) file.remove(outfile)}
  
  # now make the enclosing tabular environment
  colstring <- paste(c("l", rep("r", ncol(x))), collapse="")
  if(open_tabular == TRUE) cat(paste("\\begin{tabular}{", colstring, "}", sep=""), "\n", file = outfile, append = T)
  
  # now, haggle with getting the multicolumn header on and written to file
  cat(paste("  &  \\multicolumn{", ncol(x), "}{c}{\\underline{", rightcolshead, "}} \\\\\n", sep=""), file = outfile, append = T)
  
  # now, the second line in xf should start with \hline
  xf[2,1] <- paste("\\hline", xf[2,1])
  write.table(format(xf, digits=2), quote=F, sep="  &  ", col.names=F, row.names=F, eol="  \\\\  \n", file = outfile, append = T)
  
  # finally, enclose the tabular environment if need be
  if(close_tabular == TRUE) cat("\\end{tabular} \n", file = outfile, append = T)
}