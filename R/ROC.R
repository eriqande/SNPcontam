#' Receiver Operator Characteristic (ROC) curve for Statistical Test
#' 
#' Creates a ROC curve for a statistical test. The function can be used to create an ROC
#' for a test using heterozygosity proportion and a likelihood ratio test for contaminated
#' samples.
#' @param n An integer. Number of observations for the x variable (proportion of false positives).
#' The proportion of false positives runs from 0 to 1.
#' @param dstr0 A vector containing the values from the distribution representing the null
#' hypothesis.
#' @param dstr1 A vector containing the values from the other distribution, representing
#' the alternative hypothesis.
#' @return Returns a Receiver Operator Characteristic (ROC) curve for the test:
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
#' # x-values from 0 to 1 with steps of 0.1
#' ROC(10,ratio[[1]],ratio_c[[1]])
ROC <- function(n, dtr0, dtr1){
  fp <- seq(0,1,1/n)
  tp <- rep(0,n)
  for (i in 1:n+1){
    thresh <- threshold(dtr0, dtr1, fp[i])
    below <- dtr1[dtr1 < thresh$threshold]
    tp[i] <- 1 - length(below)/length(dtr1)
}
plot(fp,tp)
}