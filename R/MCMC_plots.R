#' @export
MCMC_plots <- function(dataframe){
# Four boxplots for Posterior means of z  
  loci_labeller <- function(var, value){
    value <- as.character(value)
    if (var=="loci_number") { 
      value[value=="20"] <- "20 Loci"
      value[value=="60"]   <- "60 Loci"
      value[value=="100"] <- "100 Loci"
      value[value=="200"] <- "200 Loci"
    }
    return(value)
  }
  p1 <- ggplot(subset(z2,contam_prob==0.025), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) +ggtitle("5 contaminated Samples") + ylab("Posterior Mean of z") + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
  p2 <- ggplot(subset(z2,contam_prob==0.075), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("15 contaminated Samples") + theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())
  p3 <- ggplot(subset(z2,contam_prob==0.2), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("40 contaminated Samples") + ylab("Posterior Mean of z") + xlab("Contaminated Status")
  p4 <- ggplot(subset(z2,contam_prob==0.5), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("100 contaminated Samples") + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) + xlab("Contaminated Status")
  
  grid.arrange(p1, p2, p3, p4, ncol=2) 
  }

# Rho Plots

# Allele Plot