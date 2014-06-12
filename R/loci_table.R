#' @export
loci_table <- function(sample_data, N, n, p, l, alpha, beta, lambda, inters, threshold){
  info <- data.frame(Loci=numeric(n), Mean_p=numeric(n),Mean_error=numeric(n),False_Positives=numeric(n),False_Negatives=numeric(n))
  for(i in 1:n){
    L <- 10*i
    MCMC <- test_MCMC(sample_data, N, L, p, l, alpha, beta, lambda, inters, threshold)
    info$Loci[i] <- L
    info$Mean_p[i] <- MCMC$mean_p
    info$Mean_error[i] <- MCMC$mean_error
    info$False_Positives[i] <- MCMC$false_pos
    info$False_Negatives[i] <- MCMC$false_neg
  }
  print(info, floating=FALSE)
}