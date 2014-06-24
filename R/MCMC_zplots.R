#' @export
MCMC_zplots <- function(z_df){
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
  # Four boxplots for Posterior means of z  
  p1 <- ggplot(subset(z_df,contam_prob==0.025), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("5 Contaminated Samples") + ylab("Posterior Mean of z") + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
  p2 <- ggplot(subset(z_df,contam_prob==0.075), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("15 Contaminated Samples") + ylab("Posterior Mean of z") + xlab("Contaminated Status")
  p3 <- ggplot(subset(z_df,contam_prob==0.2), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("40 Contaminated Samples") + ylab("Posterior Mean of z") + xlab("Contaminated Status")
  p4 <- ggplot(subset(z_df,contam_prob==0.5), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("100 Contaminated Samples") + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) + xlab("Contaminated Status")
  p5 <- ggplot(subset(z_df,contam_prob==0.0), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("0 Contaminated Samples") + ylab("Posterior MEan of z") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  a <- grid.arrange(arrangeGrob(p5, p1,nrow=1), arrangeGrob(p2, p3, p4, nrow = 1, ncol = 3), nrow=2) 
  print(a)
}