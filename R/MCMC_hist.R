# Makes 20 histograms of the z values and makes 5 graphs (4 histograms each)
# stores them in pdf documents
MCMC_hist <- function(z_df,types, width, height){
  hp <- list()
  for(i in 1:length(types)){
    p1 <- ggplot(subset(z_df,contam_prob==types[[i]]$rho & loci_number==types[[i]]$numL), aes(z,fill=as.factor(z_id))) + 
      geom_histogram(binwidth=.05, alpha = .8, position="identity") + theme_bw() +
      scale_fill_manual(values=c("grey20","grey60"), name="Contamination Status") +  
      xlab("Posterior Mean") + ggtitle(paste(toString(types[[i]]$numL),"Loci"))
    hp[[i]] <- p1
  }
  hps <- lapply(1:length(types),function(x) {hp[[x]] + theme(axis.title=element_blank(),legend.position="none")})
  
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  legend1 <- g_legend(hp[[1]])
  legend2 <- g_legend(hp[[5]])
  lwidth1 <- sum(legend1$width)
  lwidth2 <- sum(legend2$width)
  
  dev.off()
  pdf("histogram0.pdf",width=width,height=height)
  grid.arrange(arrangeGrob(hps[[1]],hps[[2]],hps[[3]],hps[[4]], left="Count"),
                   legend1,widths=unit.c(unit(1,"npc")- lwidth1, lwidth1), nrow=1)

  dev.off()
  pdf("histogram25.pdf",width=width,height=height)
  b <- grid.arrange(arrangeGrob(hps[[5]],hps[[6]],hps[[7]],hps[[8]], left="Count"),
                   legend2,widths=unit.c(unit(1,"npc")- lwidth2, lwidth2), nrow=1)
  dev.off()
  pdf("histogram75.pdf",width=width,height=height)
  c <- grid.arrange(arrangeGrob(hps[[9]],hps[[10]],hps[[11]],hps[[12]], left="Count"),
                   legend2,widths=unit.c(unit(1,"npc")- lwidth2, lwidth2), nrow=1)
  dev.off()
  pdf("histogram2.pdf",width=width,height=height)
  d <- grid.arrange(arrangeGrob(hps[[13]],hps[[14]],hps[[15]],hps[[16]], left="Count"),
                   legend2,widths=unit.c(unit(1,"npc")- lwidth2, lwidth2), nrow=1)
  dev.off()
  pdf("histogram5.pdf",width=width,height=height)
  f <- grid.arrange(arrangeGrob(hps[[17]],hps[[18]],hps[[19]],hps[[20]], left="Count"),
                   legend2,widths=unit.c(unit(1,"npc")- lwidth2, lwidth2), nrow=1)
  dev.off()
}