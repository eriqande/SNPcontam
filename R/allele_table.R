#' @export
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