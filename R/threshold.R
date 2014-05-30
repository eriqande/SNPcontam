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
  # number of non-contaminated samples that can be indentified as contaminted to obtain desired p-value
  x <- ceiling(length(dtr0)*p) 
  # dtr_a is vector with only non-contaminated values greater than or equal to the minimum contaminated value
  dtr_a <- dtr0[dtr0 >= min(dtr1)]
  # finds number of non-contaminated values greater than or equal to min. contaminated value
  lr <- length(dtr_a)
  # if lr is less than equal to x, threshold value is the min. contaminated value
  if (lr <= x){
    thresh <- min(dtr1) # sets threshold to be min. contaminated value
  } else { 
    dtr_a <- sort(dtr_a) # dtr_a sorted from low to high
    # xn used as index to find the threshold. There are x values >= the xn value of dtra
    xn <- lr - x + 1
    t <- dtr_a[xn] # possible threshold unless there are multiple values equal to t
    # x2 is vector of all dtr_a values <= t with the first value removed
    x2 <- head(dtr_a[dtr_a <= t], -1)
    # if loop checks if t appears multiple times in the non contaminated data set
    # if t appears below than it cannot be threshold because there will be more the x false positives
    if (t %in% x2){
      dtr_new <- dtr_a[dtr_a < t] # removes values of t from dtr_a
      thresh <- tail(dtr_new, n=1) # threshold is set to be largest less than t
     }else{
      thresh <- t # if t is not below than threshold is t
      }
    }
  list(threshold = thresh) # output list with threshold
}