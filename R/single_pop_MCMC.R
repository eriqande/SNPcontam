#' @export
single_pop_MCMC <- function(data,locstart = 5, inters = 1000,alpha=0.5,beta=0.5,lambda=0.5, burnin = 100){
  change_to_ones_zeros_twos <- function(data, locstart){
    locs <- colnames(data)[locstart:ncol(data)]
  
    S <- data[, locs]  # grab just the genetic data
  
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
  
    tmp <- lapply(seq(1, ncol(S2), 2), function(i) {
      rowSums(
        cbind(
          as.integer(S2[, i]) - 1,   # the minus ones here make the alleles 0 and 1, rather than 1 and 2
          as.integer(S2[, i+1]) - 1
        )
      )}
    )
  
    prepared_data <- do.call(cbind, args = tmp)
    rownames(prepared_data) <- rownames(S2)
    colnames(prepared_data) <- colnames(S2)[seq(1, ncol(S2), 2)]
    prepared_data
  }
  
  prepared_data <- t(change_to_ones_zeros_twos(data = data,locstart = locstart))
  
  MCMC <- contam_MCMC(data = prepared_data,inters = inters,alpha = alpha, beta = beta,lambda = lambda)
  
  analyze <- analyze_MCMC(MCMC = MCMC, burnin = burnin)
  
  list(MCMC_data = MCMC, analysis = analyze)
  
} 