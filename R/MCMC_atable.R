MCMC_atable <- function(allele_df,rhovals,Lvals){
  get_astats <- function(y,z){
    sub <- allele_df[allele_df$contam_prob == y & allele_df$loci_number == z,]
    a <- sub[(sub$estimates >= sub$bottomint & sub$estimates <= sub$topint),]
    total <- dim(sub)[1]
    frac <- dim(a)[1]/total
    return(frac)
  }
  values <- lapply(types,function(x) {get_astats(y = x$rho, z = x$numL)})
  mat <- round(matrix(unlist(values), nrow=length(rhovals)),4)
  rownames(mat) <- paste0(rhovals)
  colnames(mat) <- paste0(Lvals)
  return(mat)
}