#' @export
MCMC_rhoplot <- function(rho_df,rhovals){
a <- ggplot(rho_df, aes(x=factor(loci_number), y=rho_pm)) + geom_boxplot() + facet_grid(.~contam_prob) + ggtitle("Posterior Mean of Contamination Proportion") + ylab("Posterior Mean") + xlab("Loci Number")
a <- a + geom_abline(aes(intercept=rhovals,slope=0),lty="dotted")
print(a)
}