#' Full Conditional Distribution of Contamination Indicator
#' 
#' Describtion
#' @param genos 
#' @param theta 
#' @param ro
#' 
#' @return Returns a list of one named components:
#' \describe{
#'  \item{prob}{}
#' }
#' @export
#' @examples
#'
#' @export
full_z <- function(genos,theta,ro){
  like <- likelihood(theta,genos) # use likelihood function to get probability of genotype given contamination and clean
  clean <- apply(like$clean,2,prod,na.rm=TRUE)*ro # probability of non contaminated sample
  contam <- apply(like$contam,2,prod,na.rm=TRUE)*(1-ro) # probability of contaminated samples
  # full conditional probability distribution of z indicator value of 1 (indicating contamination)
  # normalized by (clean + contam) so that total probability is 1
  p <- contam/(clean + contam) 
  list(prob = p)
}
  