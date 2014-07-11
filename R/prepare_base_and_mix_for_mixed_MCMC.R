

#' Prepares a baseline and mixture file of genotypes for mixed MCMC analysis
#' 
#' This takes two data frames that should be in the same format: two columns
#' for every locus, the first column of each locus bearing a name of the locus and
#' the second one being whatever it is (but consistent between the mixture and 
#' baseline samples). Missing loci should be NA in each column of the locus.
#' Note that the columns of locus data should be the last thing in the data frame, 
#' i.e. no columns that are not genetic data should be to the right of any columns that are
#' not genetic data!  Both B and M should have rownames which are the IDs of each
#' individual. The function does this:
#' \enumerate{
#'  \item Extract loci in the mixture file from the baseline file.  So, names of
#'  loci must be consistent between the files.  Error occurs if loci in mixture file
#'  do not appear in baseline file.
#'  \item Define alleles as 0's and 1's and return counts of 0s and 1s in each baseline 
#'  population, and also return the mixture as an L x n array of 0s, 1s, and 2s.
#' }
#' @param B  Baseline data frame
#' @param B_locstart The index of the first column of genetic data in the Baseline.
#' @param B_pops  a factor vector giving the population of origin of each individual in the Baseline.  This should have levels 
#' which are ordered the way they should be.
#' @param M Mixture data frame
#' @param M_locstart  The index of the first column of genetic data in the Mixture.
#' @return This returns a list.  Letting P be the number of populations in the baseline, N be the number of individuals
#' in the mixture and L be the number of SNP loci that have exactly two alleles, this list contains the following components:
#' \describe{
#'  \item{zeros}{An L x P matrix giving the number of "0" alleles at each locus in each population.  The order of loci is as given in the 
#'  data frame M, and the order of populations is determined by the factor B_pop.}
#'  \item{one}{Same as above but for the "1" alleles.}
#'  \item{mixmat}{An L x N matrix of 0s, 1s, 2s, or NAs, giving the genotypes of the fish in the mixture.}
#' }
#' @export
#' @examples
#' # first make baseline and mixture samples
#' set.seed(5)
#' grab <- sample(1:nrow(swfsc_chinook_baseline), 400)  # grab these as a mixture
#' Base <- swfsc_chinook_baseline[-grab, ]
#' Mix <- swfsc_chinook_baseline[grab, ]
#'
#' # then prep em
#' prepped <- prepare_base_and_mix_for_mixed_MCMC(B = Base, B_locstart = 5, B_pops = Base$Pop, M = Mix, M_locstart = 5)
#'
#' names(prepped)
prepare_base_and_mix_for_mixed_MCMC <- function(B, B_locstart, B_pops, M, M_locstart) {
  
  blocs <- colnames(B)[B_locstart:ncol(B)]
  mlocs <- colnames(M)[M_locstart:ncol(M)]
  
  missing <- mlocs[!(mlocs %in% blocs)]
  if(length(missing) > 0) stop(paste("Missing these loci in baseline: ", paste(missing, collapse=", ")))
  
  B2 <- B[, mlocs]  # grab just the genetic data from M out of B
  M2 <- M[, mlocs]  # grab just the genetid data out of M
  
  BM <- rbind(B2, M2)  # bung them together to figure out which alleles there are
  
  # make a list of unique alleles
  uniq_alleles <- lapply(seq(1, ncol(BM), 2), function(x) levels(factor(c(BM[,x], BM[, x+1]))))
  names(uniq_alleles) <- colnames(BM)[seq(1, ncol(BM), 2)]
  
  # drop loci that do not have two alleles
  Have2 <- sapply(uniq_alleles, function(x) length(x)==2)
  Have2rep <- rep(Have2, each=2)  # this can be used to extract the loci from M2 and B2 that have two alleles
  
  DropTheseLoci <- names(uniq_alleles)[!Have2]
  if(length(DropTheseLoci)>0) warning(paste("Dropping these loci that have > of < then 2 alleles:", paste(DropTheseLoci, collapse=", ")))
  
  B3 <- B2[, Have2rep]
  M3 <- M2[, Have2rep]
  uniq_alleles3 <- uniq_alleles[Have2]

  
  # now, make sure that the levels of B3 and M3's alleles are appropriately set
  for(i in seq(1, ncol(B3), 2)) {
    B3[, i] <- factor(B3[, i], levels = uniq_alleles[[colnames(B3)[i]]])
    B3[, i+1] <- factor(B3[, i+1], levels = uniq_alleles[[colnames(B3)[i]]])
    M3[, i] <- factor(M3[, i], levels = uniq_alleles[[colnames(B3)[i]]])
    M3[, i+1] <- factor(M3[, i+1], levels = uniq_alleles[[colnames(B3)[i]]])
  }
  
  
  # now count the number of zero and one alleles in each population at each locus in baseline B3
  alle_counts <- lapply(seq(1, ncol(B3), 2), function(y) table(c(B3[,y], B3[,y+1]), rep(B_pops,2)))
  names(alle_counts) <- colnames(B3)[seq(1, ncol(B3), 2)]
  
  # now extract the zeros and ones matrices.  These are matrices that hold the number of 
  # "0" alleles and the number of "1" alleles in each population, respectively.
  zeros <- do.call(rbind, args = lapply(alle_counts, function(x) x[1,]))
  ones <-  do.call(rbind, args = lapply(alle_counts, function(x) x[2,]))
  
  
  # OK, now we must convert the mixture genotypes to 0s,  1s, and 2s
  tmp <- lapply(seq(1, ncol(M3), 2), function(x) {
    rowSums(
      cbind(
        as.integer(M3[, i]) - 1,   # the minus ones here make the alleles 0 and 1, rather than 1 and 2
        as.integer(M3[, i+1]) - 1
        )
      )}
    )
  
  mmat <- do.call(cbind, args = tmp)
  rownames(mmat) <- rownames(M3)
  colnames(mmat) <- colnames(M3)[seq(1, ncol(M3), 2)]
  
  
  list(zeros = zeros, ones = ones, mixmat = t(mmat))
  
}


