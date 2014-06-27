MCMC_rhoplot <- function(rho_df,rhovals){
a <- ggplot(rho_df, aes(x=factor(loci_number), y=rho_pm)) + geom_boxplot(outlier.shape=1, outlier.size=1.5) +
  facet_grid(.~contam_prob) + ylab("Posterior Mean") + xlab("Number of Loci") + geom_hline(aes(yintercept=contam_prob),lty="dotdash") +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.title.x=element_text(size=10), plot.title = element_text(hjust = 0)) + ggtitle("a")
return(a)
}