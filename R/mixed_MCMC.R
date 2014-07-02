#' @export
mixed_MCMC <- function(data, contam_data, clean_data, alpha=.5, beta=.5, inters){
  
  cbs <- nrow(contam_data)
  P <- nrow(clean_data)
  N <- ncol(data)
  L <- nrow(data)
  
  rho <- rep(0,inters+1)
  rho[1] <- rbeta(1,alpha,beta)
  u <- matrix(0,N,nrow=inters*2)
  
  pii <- matrix(0,inters+1,P)
  gm0 <- rgamma(P,shape=rep(1/P,P))
  pii[1,] <- gm0/sum(gm0)
  
  z <- matrix(0,inters,N)
  
  combos <- list()
  pops <- 1:P
  k <- 0
  for (p1 in pops) for (p2 in pops){
    k <- k + 1
    combos[[k]] <- c(p1,p2)
  }
  
  for(i in 1:inters){
    c_pi <- rep(pii[i,], each=P)*rep(pii[i,], times=P)
  # update z
    probs_c <- contam_data*c_pi
    probs0 <- clean_data*pii[i,]
    pcontam <- colSums(probs_c*rho[i])
    pclean <- colSums(probs0*(1-rho[i]))
    prob <- pcontam/(pclean + pcontam)
    z[i,] <- c(runif(N) <= prob)*1
  
  # update u
    z_c <- which(z[i,] %in% 1)
    z0 <- which(z[i,] %in% 0)
    #contaminated
    if (length(z_c) != 0){
      if (length(z_c) == 1){
        uprobs_c <- probs_c[,z_c]
        pops_c <- sample((1:cbs), size=1, prob = uprobs_c)
        }
      else {
        uprobs_c <- apply(probs_c[,z_c],2,list)
        pops_c <- sapply(uprobs_c, function(x) {sample((1:cbs),size=1, prob = unlist(x))})
        }
      u[(2*i-1):(2*i), z_c] <- unlist(combos[pops_c])
    }
    #clean
    if (length(z0) != 0){
      if (length(z0) == 1){
        uprobs0 = probs0[,z0]
        pops0 <- sample((1:P),size=1, prob = uprobs0)
      }
      else {
        uprobs0 <- apply(probs0[,z0],2,list)
        pops0 <- sapply(uprobs0, function(x) {sample((1:P), size=1, prob = unlist(x))})
      }
      u[(2*i-1),z0] <- pops0  
    }  
  
  # update pi
    clean_us <- u[1,][u[2,]==0]
    utmp <- factor(clean_us, level=(1:P))
    freqs <- as.numeric(table(utmp))
    pi_alpha <- (freqs + 1/P)
    gm <- rgamma(P,shape=(pi_alpha))
    pii[i+1,] <- gm/sum(gm)
  
  # update rho
  sum_z <- sum(z[i,]) # total number of contaminated samples
  # updates rho with derived beta distribution
  p_alpha <- sum_z + alpha # alpha parameter
  p_beta <- N - sum_z + beta # beta parameter
  rho[i+1] <- rbeta(1,p_alpha,p_beta)
  }
list(prob_contam = rho, pops = u, z = z, mixing = pii)
}