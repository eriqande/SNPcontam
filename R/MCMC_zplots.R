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
  p1 <- ggplot(subset(z_df,contam_prob==0.025), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("5 contaminated Samples") + ylab("Posterior Mean of z") + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
  p2 <- ggplot(subset(z_df,contam_prob==0.075), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("15 contaminated Samples") + theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())
  p3 <- ggplot(subset(z_df,contam_prob==0.2), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("40 contaminated Samples") + ylab("Posterior Mean of z") + xlab("Contaminated Status")
  p4 <- ggplot(subset(z_df,contam_prob==0.5), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("100 contaminated Samples") + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) + xlab("Contaminated Status")
  a <- grid.arrange(p1, p2, p3, p4, ncol=2) 
  
  d1 <- ggplot(subset(z_df,contam_prob==0.025), aes(x=factor(loci_number),y=z)) + geom_boxplot(aes(fill=factor(z_id))) + ggtitle("5 contaminated Samples") + ylab("Posterior Mean of z") + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="none")
  d2 <- ggplot(subset(z_df,contam_prob==0.075), aes(x=factor(loci_number),y=z)) + geom_boxplot(aes(fill=factor(z_id))) + ggtitle("15 contaminated Samples") + theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank()) + scale_fill_discrete(name="Contamination Status",breaks=c("0", "1"), labels=c("Non-contaminated", "Contaminated"))
  d3 <- ggplot(subset(z_df,contam_prob==0.2), aes(x=factor(loci_number),y=z)) + geom_boxplot(aes(fill=factor(z_id))) + ggtitle("40 contaminated Samples") + ylab("Posterior Mean of z") + xlab("Contaminated Status") + theme(legend.position = "none")
  d4 <- ggplot(subset(z_df,contam_prob==0.5), aes(x=factor(loci_number),y=z)) + geom_boxplot(aes(fill=factor(z_id))) + ggtitle("100 contaminated Samples") + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(), legend.position = "none") + xlab("Contaminated Status")
  b <- grid.arrange(d1, d2, d3, d4, ncol=2) 
  print(a,b)
}