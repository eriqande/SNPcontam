MCMC_alleleplot <- function(allele_df,rho,loci){
  # Allele Plot
  a <- ggplot(subset(allele_df,contam_prob == rho & loci_number==loci), aes(x=alle_freq,y=estimates)) + 
    geom_point(shape=19, size=.60) + geom_abline(aes(intercept=0,slope=1), lty="dashed") + 
    xlab("True Allele Frequency") + ylab("Posterior Mean of Allele Frequency") + 
    theme(text = element_text(size=10), axis.text=element_text(size=8), plot.title = element_text(hjust = 0)) + ggtitle("b")
  return(a)
}
# smaller titles