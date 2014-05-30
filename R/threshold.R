#' Threshold Value for False Positives
#' 
#' Computes the threshold value for a text given two different distributions (i.e. a
#' distribution of likelihood ratios from a contaminated population and a non 
#' contaminated population)
#' @param dstr0 A vector containing the values from the distribution representing the null
#' hypothesis.
#' @param dstr1 A vector containing the values from the other distribution, representing
#' the alternative hypothesis.
#' @param p A number the represents the proportion of false positives desired. For example,
#' if the dstr0 sample contains 100 elements and only 2 false positives are allowed, p = 0.2.
#' @return Returns a list of one named component:
#' \describe{
#'  \item{threshold}{A number that represents the lower bound of the rejection region of the
#'  the null hypothesis in favor of the alternative. With regards to contamination, all samples
#'  contain a test statistic greater than this number will be considered contaminated}
#' }
#' @export
#' @examples
#' # calling after running ratio tests for clean and contaminated samples
#' af <- runif(10,0,1)
#' Ran <- random_gene(10,af)
#' L <- likelihood(af,Ran$rclean)
#' LC <- likelihood(af, Ran$rcontam)
#' ratio <- lratio(L$clean, L$contam,af)
#' ratio_c <- lratio(LC$clean, LC$contam,af)
#' # non contamined ratios are null hypothesis and contaminated are alternative
#' # proportion of false positives is 0.01
#' threshold(ratio[[1]],ratio_c[[1]],0.01)
threshold <- function(dtr0, dtr1, p){
  x <- ceiling(length(dtr0)*p)
  dtr_u <- dtr0[dtr0 >= min(dtr1)]
  lr <- length(dtr_u)
  if (lr <= x){
    thresh <- min(dtr1)
  } else { 
    xn <- length(dtr_u) - x
    dtr_u <- sort(dtr_u)
    t <- dtr_u[xn]
    x2 <- head(dtr_u[dtr_u <= t], -1)
    if (length(dtr_u) != 1 & t %in% x2){
      dtr_new <- dtr_u[dtr_u < t]
      t2 <- tail(dtr_new, n=1)
      thresh <- t2
    } else{
      thresh <- t
      }
    }
  list(threshold = thresh)
}