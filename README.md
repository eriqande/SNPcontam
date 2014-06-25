# SNPcontam

This is an R package with some code implementing methods under development
to detect contaminated samples from SNP genotyes sampled in molecular
ecology contexts.  

Currently, the directory `not-package` contains things that will not go into
the final R package.  This includes the manuscript that we will be preparing
and several readings and resources.  

To typeset the manuscript (which is really just a stub at this point)
`cd` into the `not-package/manuscript` directory and issue this
command (note that you have to have a TeX system installed on your
system):
```
latexmk -pdf snp-contam-main
```
The resulting PDF file, `snp-contam-main` should then be complete.

## Preparing Data for the Markov Chain Monte Carlo (MCMC) Method
The MCMC function requires that the input data be a matrix with columns equal to number individuals and rows equal to number of loci where each element is the genotype of an individual at a particular locus.  Genotypes are given by the number of "1" alleles and should be 0s, 1s, 2s, or NA.  Functions from the package `fullsniplings` can be used to turn a data frame of SNP alleles into values of 0, 1, 2, or NA using the code below:
```
library(fullsniplings)
snp_genos <- get_snp_genos(tmp_data)
data <- snp_genos$mat
```
Any columns or rows contain non-SNP information (i.e. information regarding the origin of the individual `tmp_data` should only be removed before using the function `get_snp_genos`.

## Running the Code for the Markov Chain Monte Carlo (MCMC) Method

The functions below are used to run the MCMC for a set of data called `data`.
```
MCMC <- contam_MCMC(data)
means <- analyze_MCMC(MCMC)
```
The function`contam_MCMC` performs the MCMC method on the data, returning the all the values of the allele frequencies, proportion of contaminated samples, and z values for each sweep. If the input into `contam_MCMC` is only `data`, the function will use default values of `inters = 1000`, `alpha = 0.5`, `beta=0.5`, and `lambda=0.5`.

The function `analyze_MCMC` computes the posterior means of the allele frequencies, contamination proportion, and z values.  If the MCMC results are the only input data, `analyze_MCMC` will use a default value of `burnin=100`.

## Terms 

As a work partially of the United States Government, this package is in the
public domain within the United States. 

Of course, it would be awfully lame of anyone to take anything found within
here and use it as their own before we tried publishing any of this, should
we choose to do that.

See TERMS.md for more information.

