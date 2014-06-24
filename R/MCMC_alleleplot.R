#' @export
MCMC_alleleplot <- function(allele_df,rho,loci){
  # Allele Plot
  a <- ggplot(subset(allele_df,contam_prob == rho & loci_number==loci), aes(x=alle_freq,y=estimates)) + geom_point(shape=19) + geom_abline(aes(intercept=0,slope=1), lty="dashed") + xlab("True Allele Frequency") + ylab("Posterior Mean of Allele Frequencies") + ggtitle("Allele Frequency Comparision") + geom_errorbar(width=.0075, aes(ymin=bottomint, ymax=topint))
  print(a)
}