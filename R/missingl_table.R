missingl_table <- function(sample_data, N, L, p, n, alpha, beta, lambda, inters, threshold){
  info <- data.frame(Contam_Loci=numeric(n), Mean_p=numeric(n),Mean_error=numeric(n),False_Positives=numeric(n),False_Negatives=numeric(n))
  for(i in 1:n){
    l <- 1.1 - .1*i
    MCMC <- test_MCMC(sample_data, N, L, p, l, alpha, beta, lambda, inters, threshold)
    info$Contam_Loci[i] <- l
    info$Mean_p[i] <- MCMC$mean_p
    info$Mean_error[i] <- MCMC$mean_error
    info$False_Positives[i] <- MCMC$false_pos
    info$False_Negatives[i] <- MCMC$false_neg
  }
  print(info, floating=FALSE)
}