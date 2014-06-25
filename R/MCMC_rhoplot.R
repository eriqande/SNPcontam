#' @export
MCMC_rhoplot <- function(rho_df,rhovals){
rho_df$rhov <- rhovals
a <- ggplot(rho_df, aes(x=factor(loci_number), y=rho_pm)) + geom_boxplot() + facet_grid(.~contam_prob) + ggtitle("Posterior Mean of Contamination Proportion") + ylab("Posterior Mean") + xlab("Loci Number") + geom_hline(aes(yintercept=contam_prob),lty="dotdash")
print(a)
}