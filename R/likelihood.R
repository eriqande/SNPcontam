#Likelihood of genotype given contamination and noncontamination
likelihood <- function(gt,af) {
  #calculate probablitiy for noncontaminated
  g0 <- af^2; g1 <- 2*af*(1-af); g2 <- (1-af)^2
  #calculate probability for contaminated
  gc0 <- af^4; gc1 <- 1 - (af^4 + (1-af)^4); gc2 <- (1-af)^4 
  #make vector of likelihood for each non contaminated gene
  a <- gt==0; b <- gt==1; c <- gt==2
  g <- g0*a + g1*b + g2*c
  #make vector of likelihood for each contaminated gene
  ac <- gt==0; bc <- gt==1; cc <- gt==2
  gc <- gc0*ac + gc1*bc + gc2*cc
  list(contam = gc, clean = g)
}
