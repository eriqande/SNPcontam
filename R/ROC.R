#' Receiver Operator Characteristic (ROC) curve for Statistical Test
#' 
#' Creates an ROC curve for a statistical test. With regards to contamination, the function 
#' is used to create an ROC for a test using heterozygosity proportion and a likelihood ratio test for contaminated
#' samples.
#' @param n An integer that is the number of observations for the x variable (proportion of 
#' false positives). The proportion of false positives runs from 0 to 1, and 1/n determines
#' the step between each x-value.
#' @param dstr0 A vector containing the values from the distribution representing the null
#' hypothesis, which in this case is the distribution created from the non-contaminated genotypes.
#' @param dstr1 A vector containing the values from the other distribution, representing
#' the alternative hypothesis. In this case, dstr1 is contains the values for the distribution
#' created using contamnited genotypes.
#' @param title A string that will be the main title of the ROC graph.
#' @return Returns a Receiver Operator Characteristic (ROC) curve for the test:
#' \describe{
#'  The x-axis is the proportion of false positives.  The y-axis is the proportion of
#'  true positives (1 - false negatives).  the x- and y-axis run from 0 to 1.
#' }
#' @export
#' @examples
#' # call after running ratio tests
#' af <- runif(10,0,1)
#' Ran <- random_gene(10,af)
#' L <- likelihood(af,Ran$rclean)
#' LC <- likelihood(af, Ran$rcontam)
#' ratio <- lratio(L$clean, L$contam,af)
#' ratio_c <- lratio(LC$clean, LC$contam,af)
#' # x-values from 0 to 1 with steps of 0.1
#' ROC(10,ratio[[1]],ratio_c[[1]],"ROC Curve for Likelihood Ratio Test")
#' 
#' #call after running hetero function for ROC curve of heterozygosity test
#' af <- runif(10,0,1)
#' Ran <- random_gene(10,af)
#' hetero <- hetero(Ran$rclean,af)
#' hetero_c <- hetero(Ran$rcontam, af)
#' #x-values from 0 to 1 with steps of 0.1
#' ROC(10,hetero[[1]],hetero_c[[1]],"ROC Curve for Heterozygosity Proprotion Test")
ROC <- function(n, dtr0, dtr1,title){
  fp <- seq(0,1,1/n) # values of false positives proportion
  tp <- rep(0,n) # empty vector for true positives to be recorded
  for (i in 1:n+1){
    # calss threshold fucntion to determine threshold for false positives
    thresh <- threshold(dtr0, dtr1, fp[i])
    # vector containing all values less than the threshold
    # these are false negatives because they are contaminated but test says they are not
    below <- dtr1[dtr1 < thresh$threshold]
    # true positives proportion are equal to 1 - false negatives proportion
    tp[i] <- 1 - length(below)/length(dtr1)
}
plot(fp,tp,xlab = "Proportion of False Positives", ylab = "Proportion of True Positives", main = title) # plots true positives vs. false positives
}