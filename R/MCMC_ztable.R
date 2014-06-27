MCMC_ztable <- function(z_df,types,rhovals,Lvals){
  get_zstats <- function(y,z){
    sub <- z_df[z_df$contam_prob==y & z_df$loci_number==z,]
    t1 <- length(sub$z[sub$z_id==1])
    t0 <- length(sub$z[sub$z >= .9])
    frac_1 <- length(sub$z[sub$z_id==1 & sub$z >= .9])/t1
    frac_0 <- length(sub$z[sub$z_id==0 & sub$z >= .9])/t0
    list(frac_1,frac_0)
  }
  tmp <- lapply(types,function(x) {get_zstats(y = x$rho, z = x$numL)})
  values <- lapply(tmp, function(x) {if (is.na(x[[1]])) {x[[1]] = NA}; if (is.na(x[[2]])){x[[2]] = 0}; return(x)})
  mat <- round(matrix(unlist(values), nrow=length(rhovals)),4)
  rownames(mat) <- paste0(rhovals)
  colnames(mat) <- paste0(rep(Lvals,each=2)," Loci")
  return(mat)
}