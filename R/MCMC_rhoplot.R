#' @export
MCMC_rhoplot <- function(rho_df){
a <- ggplot(rho_df, aes(x=factor(loci_number), y=rho_pm)) + geom_boxplot() + facet_grid(.~contam_prob) + ggtitle("Posterior Mean of Contamination Proportion") + ylab("Posterior Mean") + xlab("Loci Number") + geom_abline(aes(intercept=0.025,slope=0),lty="dotted") + geom_abline(aes(intercept=0.075, slope=0), lty="dotted") + geom_abline(aes(intercept=0.2,slope=0),lty="dotted") + geom_abline(aes(intercept=0.5,slope=0), lty="dotted") + scale_y_continuous(breaks = seq(0, 0.55, 0.025))
print(a)
}