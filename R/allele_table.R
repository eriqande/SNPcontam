#' @export
allele_table <- function(data,Lvals,rhovals){
  nL <- length(Lvals)
  nrho <- length(rhovals)
  data2 <- data[,-(5:6)]
  data2$difference <- abs(data2$alle_freq - data2$estimates)
  rho <- numeric(0)
  Ls <- numeric(0)
  for (r in rhovals){
    a1 <- rep(r,nL)
    rho <- c(rho,a1)
  }
  for (L in Lvals){
    a2 <- L
    Ls <- c(Ls,a2)
  }
  table_data <- data.frame(Rho_Value = rho, Number_of_Loci= rep(Ls,nrho), Average_Difference = rep(0,nL*nrho),)
  for (i in 1:(nL*nrho)){
    loci <- table_data$Number_of_Loci[i]
    r <- table_data$Rho_Value[i]
    table_data$Average_Difference[i] = mean(data2$difference[data2$loci_number == loci & data2$rho_value == r])
  }
  colnames(table_data) <- c("Proportion of Contaminated Samples", "Number of Loci", "Mean Absolute Difference")
  print(xtable(table_data,digits=c(0,3,0,3)),include.rownames=FALSE)
}