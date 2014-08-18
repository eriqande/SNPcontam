# Functions Needed for Running the Simulation
## Simulate Random Genotypes
simulate_genos <- function(N,L,p,l, afreqs){
  n_contam <- round(N*p,0)
  l_contam <- round(L*l,0)
  g <- matrix(0,L,N)
  # get genotype frequencies
  all_L <- ncol(afreqs)
  loci <- sample(1:all_L,L,replace = TRUE)
  genes <- afreqs[2,loci] + runif(L,-.01,.01)
  genes[genes==0.0] <- .01
  genes[genes<0.0] <- -genes[genes<0.0]
  genes[genes>1.0] <- 2.0 - genes[genes>1.0]
  gfreqs <- likelihood(genes)
  
  # randomly pick individuals to be contaminated
  if (p==0){
    for (i in 1:L){
      g[i,] <- sample(0:2,N,TRUE,prob=gfreqs$clean_prob[,i])
    }
  } else{
    id <- sample(1:N,n_contam,replace = FALSE)
    contam_loci <- sample(1:L,l_contam,replace = FALSE)
    
    #simulate genotypes
    for (i in 1:L){
      if (i %in% contam_loci){
        # stores genotype at locus in all individuals
        g[i,id] <- sample(0:2, n_contam, TRUE, prob = gfreqs$contam_prob[,i])
        g[i,-(id)] <- sample(0:2,(N-n_contam),TRUE,prob = gfreqs$clean_prob[,i])
      } else{
        g[i,] <- sample(0:2,N,TRUE,prob=gfreqs$clean_prob[,i])
      }
    }
  }
  if (p==0){
    contam_id <- 0
  } else{
    contam_id = id
  }
  list(geno = g, contam_id = contam_id, loci = loci, afreqs = genes)
}

## test_MCMC is needed ot collect some of the data on the performance of the MCMC
test_MCMC <- function(afreqs, N, L, p, l=1, alpha=0.5, beta=0.5, lambda=0.5, inters, burnin=100){
  real_p <- p
  sim <- simulate_genos(N,L,p,l,afreqs)
  MCMC <- contam_MCMC(sim$geno,inters,alpha,beta,lambda)
  
  #compare allele frequencies
  alleles <- MCMC$allele_freq[-(1:burnin),]
  a_means <- colMeans(alleles) 
  bottom <- apply(alleles, 2, quantile, probs = 0.05)
  top <- apply(alleles,2, quantile, probs = 0.95)
  
  # calculate mean proportion contaminated
  p2 = mean(MCMC$prob_contam[-(1:burnin)])
  
  # compare contaminated individuals
  z <- MCMC$z[-(1:burnin),]
  z_pm <- colMeans(z)
  z_id <- rep(0,N)
  z_id[sim$contam_id] = 1
  
  list(afreqs = sim$afreqs, bottom_int = bottom, top_int = top , alle_pm = a_means ,z_pm = z_pm , z_id = z_id, rho_pm = p2 )
  
}

## Running the Simulation
MCMC_sims <- function(sample_data, N, Lvals, rhovals, l=1, alpha=0.5, beta=0.5, lambda=0.5, inters=1000,n){
  get_allele_freqs <- function(data, locstart=5){
    locs <- colnames(data)[locstart:ncol(data)]
    
    S <- data[, locs]  # grab just the genetic data
    
    uniq_alleles <- lapply(seq(1, ncol(S), 2), function(x) levels(factor(c(S[,x], S[, x+1]))))
    names(uniq_alleles) <- colnames(S)[seq(1, ncol(S), 2)]
    
    # drop loci that do not have two alleles
    Have2 <- sapply(uniq_alleles, function(x) length(x)==2)
    Have2rep <- rep(Have2, each=2)  # this can be used to extract the loci from M2 and B2 that have two alleles
    DropTheseLoci <- names(uniq_alleles)[!Have2]
    if(length(DropTheseLoci)>0) warning(paste("Dropping these loci that have > of < then 2 alleles:", paste(DropTheseLoci, collapse=", ")))
    
    S2 <- S[, Have2rep]
    uniq_alleles2 <- uniq_alleles[Have2]
    
    # now, make sure that the levels of S2's alleles are appropriately set
    for(i in seq(1, ncol(S2), 2)) {
      S2[, i] <- factor(S2[, i], levels = uniq_alleles2[[colnames(S2)[i]]])
      S2[, i+1] <- factor(S2[, i+1], levels = uniq_alleles2[[colnames(S2)[i]]])
    }
    #get zeros and one allele counts
    alle_counts <- lapply(seq(1, ncol(S2), 2), function(y) table(c(S2[, y], S2[, y+1])))
    names(alle_counts) <- colnames(S2)[seq(1, ncol(S2), 2)]
    
    # now extract the zeros and ones matrices.  These are matrices that hold the number of 
    # "0" alleles and the number of "1" alleles in each population, respectively.
    zeros <- do.call(cbind, args = lapply(alle_counts, function(x) x[1]/sum(x)))
    ones <-  do.call(cbind, args = lapply(alle_counts, function(x) x[2]/sum(x)))
    allele_freqs <- rbind(zeros,ones)
  }
  afreqs <- get_allele_freqs(sample_data)
  
  
  
  # create list with all simulations
  types <- list() 
  k <- 0
  num_Ls <- length(Lvals) # number of Loci values
  num_rho <- length(rhovals) # number of rho values
  # for loop that sets up the different sets of trials (all combinations of Lvals and rhovals)
  for (r in rhovals) for (L in Lvals){
    k <- k+1 
    types[[k]] <- list() 
    types[[k]]$rho = r 
    types[[k]]$numL = L
  }
  
  # Runs all MCMCs.  
  # MCMC contains the z posterior means, the rho posterior means, the allele frequencies, the allele estimates, and the upper and lower bound for the 90% allele intervals
  MCMC <- mclapply(types, function(x) {
    lapply(1:n, function(y) {
      list(params = x, output = test_MCMC(afreqs = afreqs, N = N, l = l, L = x$numL, p = x$rho, inters = inters))
    })
  },mc.cores=length(types))
  
  # Makes the z output data frame
  slurp_mcmc_z_output <- function(y) {
    data.frame(z = y$output$z_pm, z_id = y$output$z_id, contam_prob = y$params$rho, loci_number = y$params$numL )
  }
  ztmp1 <- mclapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x) {df <- slurp_mcmc_z_output(x); df$rep_num=rep; df }) }) 
  ztmp2 <- unlist(ztmp1, recursive = FALSE)
  z_df <- do.call(what = rbind, args = ztmp2)
  
  # Makes the rho output data frame
  slurp_mcmc_rho_output <- function(y){
    data.frame(rho_pm = y$output$rho_pm, contam_prob = y$params$rho, loci_number = y$params$numL)
  }
  rtmp1 <- mclapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x) {df <- slurp_mcmc_rho_output(x); df$rep_num=rep; df})})
  rtmp2 <- unlist(rtmp1, recursive = FALSE)
  rho_df <- do.call(what = rbind, args = rtmp2)
  
  # Makes the allele output data frame
  slurp_mcmc_allele_output <- function(y){
    data.frame(alle_freq = y$output$afreqs, estimates = y$output$alle_pm, topint = y$output$top_int, bottomint = y$output$bottom_int, contam_prob = y$params$rho, loci_number = y$params$numL)
  }
  atmp1 <- mclapply(1:length(MCMC), function(rep) {lapply(MCMC[[rep]], function(x){df <- slurp_mcmc_allele_output(x); df$rep_num=rep; df})})
  atmp2 <- unlist(atmp1, recursive = FALSE)
  allele_df <- do.call(what = rbind, args = atmp2)
  
  list(types = types, z = z_df, allele = allele_df, rho = rho_df, rhovals = rhovals, Lvals = Lvals)
}

# Functions Needed for Plotting
## Creates the allele plot
MCMC_alleleplot <- function(allele_df,rho,loci){
  # Allele Plot
  a <- ggplot(subset(allele_df,contam_prob == rho & loci_number==loci), aes(x=alle_freq,y=estimates)) + 
    geom_point(shape=19, size=.60) + geom_abline(aes(intercept=0,slope=1), lty="dashed") + 
    xlab("True Allele Frequency") + ylab("Posterior Mean of Allele Frequency") + 
    theme(text = element_text(size=10), axis.text=element_text(size=8), plot.title = element_text(hjust = 0)) + ggtitle("b")
  return(a)
}

## Creates a boxplot for the z values
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
  p1 <- ggplot(subset(z_df,contam_prob==0.025), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("5 Contaminated Samples") + ylab("Posterior Mean of z") + theme(axis.title.x=element_blank())
  p2 <- ggplot(subset(z_df,contam_prob==0.075), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("15 Contaminated Samples") + ylab("Posterior Mean of z") + xlab("Contaminated Status")
  p3 <- ggplot(subset(z_df,contam_prob==0.2), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("40 Contaminated Samples") + xlab("Contaminated Status") + theme(axis.title.y=element_blank())
  p4 <- ggplot(subset(z_df,contam_prob==0.5), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("100 Contaminated Samples") + theme(axis.title.y=element_blank()) + xlab("Contaminated Status")
  p5 <- ggplot(subset(z_df,contam_prob==0.0), aes(x=factor(z_id),y=z)) + geom_boxplot() + facet_grid(.~loci_number, labeller=loci_labeller) + ggtitle("0 Contaminated Samples") + ylab("Posterior MEan of z") + theme(axis.title.x=element_blank())
  a <- grid.arrange(arrangeGrob(p5, p1,nrow=1), arrangeGrob(p2, p3, p4, nrow = 1, ncol = 3), nrow=2)
  print(a)
}

## Creates a boxplots for the rho values
MCMC_rhoplot <- function(rho_df,rhovals){
  a <- ggplot(rho_df, aes(x=factor(loci_number), y=rho_pm)) + geom_boxplot(outlier.shape=1, outlier.size=1.5) +
    facet_grid(.~contam_prob) + ylab("Posterior Mean") + xlab("Number of Loci") + geom_hline(aes(yintercept=contam_prob),lty="dotdash") +
    theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.title.x=element_text(size=10), plot.title = element_text(hjust = 0)) + ggtitle("a")
  return(a)
}

## Creates histograms for the rho z values
MCMC_hist <- function(z_df,types, width, height, outpath){
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
  pdf(file.path(outpath, "histogram_0.pdf"), width=width,height=height)
  grid.arrange(arrangeGrob(hps[[1]],hps[[2]],hps[[3]],hps[[4]], left="Count"),
               legend1,widths=unit.c(unit(1,"npc")- lwidth1, lwidth1), nrow=1)
  
  dev.off()
  pdf(file.path(outpath, "histogram_0.025.pdf"), width=width,height=height)
  b <- grid.arrange(arrangeGrob(hps[[5]],hps[[6]],hps[[7]],hps[[8]], left="Count"),
                    legend2,widths=unit.c(unit(1,"npc")- lwidth2, lwidth2), nrow=1)
  dev.off()
  pdf(file.path(outpath, "histogram_0.075.pdf"), width=width,height=height)
  c <- grid.arrange(arrangeGrob(hps[[9]],hps[[10]],hps[[11]],hps[[12]], left="Count"),
                    legend2,widths=unit.c(unit(1,"npc")- lwidth2, lwidth2), nrow=1)
  dev.off()
  pdf(file.path(outpath, "histogram_0.2.pdf"), width=width,height=height)
  d <- grid.arrange(arrangeGrob(hps[[13]],hps[[14]],hps[[15]],hps[[16]], left="Count"),
                    legend2,widths=unit.c(unit(1,"npc")- lwidth2, lwidth2), nrow=1)
  dev.off()
  pdf(file.path(outpath, "histogram_0.5.pdf"), width=width,height=height)
  f <- grid.arrange(arrangeGrob(hps[[17]],hps[[18]],hps[[19]],hps[[20]], left="Count"),
                    legend2,widths=unit.c(unit(1,"npc")- lwidth2, lwidth2), nrow=1)
  dev.off()
}

## Creates a number ofMCMC_allplots <- function(sim,rho,loci){
MCMC_allplots <- function(sim,rho,loci){
  MCMC_zplots(sim$z) # plots z values
  MCMC_alleleplot(sim$allele,rho,loci) # plots allele data
  MCMC_rhoplot(sim$rho, sim$rhovals) # plots rho data
  allele_table(sim$types,sim$allele,sim$Lvals,sim$rhovals)
}

# Functions Needed for Making Tables
## Creates a table with absolute mean difference of allele difference
allele_table <- function(types, data, Lvals, rhovals){
  data$difference <- abs(data$alle_freq - data$estimates) # finds absolute difference between estimate and real value fo allele frequency
  # loops set up rho and L values for the data frame
  get_output <- function(y){
    rho <- y$rho
    loci <- y$numL
    df <- data.frame(Rho_Value = rho, Number_of_Loci = loci)
    df$Average_Difference <- mean(data$difference[data$loci_number == loci & data$contam_prob == rho])
    return(df)  
  }
  tmp <- lapply(types, function(x) {df <- get_output(x); df})
  table_data <- do.call(what = rbind, args = tmp)
  
  
  # get rid of extra rho values
  x <- table_data$Rho_Value
  reps <- c(FALSE,x[-1]==x[-length(x)])
  table_data$Rho_Value[reps] <- NA
  
  #change column names
  colnames(table_data) <- c("Rho Value", "Number of Loci", "Mean Absolute Difference")
  
  #print table
  nrho <- length(rhovals)
  nL <- length(Lvals)
  lines <- numeric(0)
  for (i in 1:nrho){
    lines <- c(lines,i)
  }
  tab <- xtable(table_data,digits=c(0,3,0,3))
  align(tab) <- "c|r|r|r|"
  print(tab,include.rownames=FALSE,hline.after=c(-1,0,nL*lines))
}

## Creates a matrix with information for the z value table
## z is the z_output data frame from the simulation.  PPlim
## is the posterior prob cutoff for which you call something
## contaminated.  We will use 0.5.
MCMC_ztable <- function(z, PPlim = 0.5) {
  # NP <== "No Posterior"  (We do not include posterior >0.9 as a factor)
  NP <- table(list(rho=z$contam_prob, L=z$loci_number, Z=z$z_id) )
  
  # WP <== "with posterior >0.9" included as a factor
  WP <- table(list(rho=z$contam_prob, L=z$loci_number, Z=z$z_id, GreaterThanPPlim=z$z>PPlim) )
  
  # uncontams with post prob > PPlim
  FalsePos <- WP[,,"0","TRUE"] / NP[,,"0"]
  
  # true contams with post prob > PPlim
  TruePos <- WP[,,"1","TRUE"] / NP[,,"1"]
  
  list(FP = FalsePos, TP = TruePos)
}

make_tabular <- function(x, dd=2, leftcolhead, rightcolshead, outfile="", open_tabular = TRUE, close_tabular = TRUE, append = FALSE) {
  nr <- nrow(x)
  nc <- ncol(x)
  xf <- format(round(x, digits=dd), digits=dd)  # make them strings formatted as desired
  xf <- cbind(rownames(xf), xf)  # add the rownames as a column
  xf <- rbind(c(leftcolhead, colnames(x)), xf) # add the colnames as the top row
  
  
  # remove outfile if it exists
  if(append == FALSE) {if(file.exists(outfile)) file.remove(outfile)}
  
  # now make the enclosing tabular environment
  colstring <- paste(c("l", rep("r", ncol(x))), collapse="")
  if(open_tabular == TRUE) cat(paste("\\begin{tabular}{", colstring, "}", sep=""), "\n", file = outfile, append = T)
  
  # now, haggle with getting the multicolumn header on and written to file
  cat(paste("  &  \\multicolumn{", ncol(x), "}{c}{\\underline{", rightcolshead, "}} \\\\\n", sep=""), file = outfile, append = T)
  
  # now, the second line in xf should start with \hline
  xf[2,1] <- paste("\\hline", xf[2,1])
  write.table(format(xf, digits=2), quote=F, sep="  &  ", col.names=F, row.names=F, eol="  \\\\  \n", file = outfile, append = T)
  
  # finally, enclose the tabular environment if need be
  if(close_tabular == TRUE) cat("\\end{tabular} \n", file = outfile, append = T)
}

