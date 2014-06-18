#' @export
MCMC_alleleplot <- function(allele_df){
  # Allele Plot
  a <- ggplot(subset(allele_df,rho_value == 0.2 & loci_number==20), aes(x=alle_freq,y=estimates)) + geom_point(shape=19) + geom_abline(aes(intercept=0,slope=1), lty="dashed") + xlab("True Allele Frequency") + ylab("Posterior Mean of Allele Frequencies") + ggtitle("Allele Frequency Comparision") + geom_errorbar(width=.0075, aes(ymin=bottomint, ymax=topint))
  print(a)
}