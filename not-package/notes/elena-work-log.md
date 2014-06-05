---
title: "Elena's Worklog"
author: "Elena Venable"
output:
  html_document:
    toc: true
    theme: united
  pdf_document:
    toc: true
---

## 5/27/14
### Setting up Git Repository
* Set up a Github account under the username "evenable".
* Set up a folder "Hollings" in Terminal
    * `cd` moves to home directory
    * `mkdir Hollings` makes Hollings directory
* Cloned the SNPcontam repository using the url https://github.com/eriqande/SNPcontam.git.
    * `git clone https://github.com/eriqande/SNPcontam.git` clones the SNPcontam repository to my Hollings directory
* Changed some of the git configurations
    * `git config --global user.name "evenable"` to input name
    * `git config --global user.email Elena_venable@brown.edu`: to input emai address
    * `git config --global core.editor emacs` changed the default editor to emacs
* Other important git commands
    * `git status` 
    * `git push origin master` pushes my local changes to SNPcontam repository
* Trying to get password caching to work
    * We tried this: `git config --global credential.helper osxkeychain`
    * Then we committed and pushed, and this cached the password

### Work in Rstudio
* Downloaded packages needed for SNPcontam project
    * `knitr`, `devtools`, and `roxygen2`
* Git can also be accessed and used from Rstudio undert he "Git" tab

### R code for Likelihood
* Wrote R function for determining likelihood of genotypes if contaminated or not contaminated given allele frequencies
    * the function: `likelihood(gt,af)`
    * `gt` is genotypes as 0s, 1s, and 2s, where 1s and 2s are homozygous and 1s are heterozygous
    * `af` is the allele frequencies of the loci
  
## 5/28/14
### Work Plan
* Began to formulate Work Plan
      * Likelihood Ratio Test
          * $\Lambda= log(\frac{p(\theta|con)}{p(\theta | noncon)})$ 
          * Question: How much better is a likelihood ratio test than heterozygousity test?
              * Evaluate using the Area Under the Curve (AUC) of a Receiver Operator Characteristic curve (ROC)
          * ROC Curve: True Positives (1 - False Negatives) vs. False Positives
      * Lab experiment with DNA

### Work in R
* To do list
    * Code multiplying together loci to get one number of likelihood prob. for each individual.
    * Code calculating likelihood ratio for individual
    * Code to randomly generate contaminated and non contaminated distribution
    * Code to determine heterozygousity ratio
    * Plot True + vs. False + (ROC curve)
    * Find AUC of ROC
  
## 5/29/14
### Work in R
* R Coding: Finalized Drafts of Codes
    + `ROC`: creating ROC curve
    + `threshold`: finding threshold likelihood ratio/heterozygousity ratio for p-value
    + `lratio`: calculating likelihood ratio genotypes
    + `hetero`: calculating proportion of heterozygous loci
    + `random_gene`: creates n random genotypes using genotype frequencies
* R Documention
    + Wrote the documentation for `lratio`, `random_gene`, and `hetero`
    
## 5/30/14
### Work in R
* Documentation of `ROC` and `threshold` although not sure if they work properly
* Wrote code to run all functions with random genes to build ROC curves and histograms
```
af <- runif(20,0,1)
Ran <- random_gene(1000,af)
L <- likelihood(af,Ran$rclean)
LC <- likelihood(af, Ran$rcontam)
ratio <- lratio(L$clean, L$contam,af)
ratio_c <- lratio(LC$clean, LC$contam,af)
hetero <- hetero(Ran$rclean,af)
hetero_c <- hetero(Ran$rcontam, af)
par(mfrow=c(2,2))
ROC(100,ratio[[1]], ratio_c[[1]])
hist(ratio[[1]],col=rgb(1,0,0,0.5),xlim = c(-40,20),main="Likelihood Ratio Histogram")
hist(ratio_c[[1]],col=rgb(0,0,1,0.5),add=T)
ROC(100,hetero[[1]],hetero_c[[1]])
hist(hetero[[1]],col=rgb(1,0,0,0.5),xlim = c(0,1),main="Heterozygousity Histogram")
hist(hetero_c[[1]],col=rgb(0,0,1,0.5),add=T)
```
* Appeared as if likelihood not better than heterozygousity, which makes me question if codes are working correctly
* Also above 45 randomly generated alleles, the likelihood ratios (and heterozygousity proportions as well) of the contaminated and non contaminated do not overlap.
* Committed code

### Plans for Next Week
* Document and comment code
* Think about how likelihood ratio mathematically relates to number of heterozygous loci
* Literature Review of Genotyping/contaminated samples
    * Search for info on how others have dealt with contamination problems and identifiying contamination, etc.
* Begin learning about MCMC methods
    * Eric's notes can be found here : http://users.soe.ucsc.edu/~eriq/dokuwiki/doku.php?id=sisg:sisg
* Think about how I would like to write up the results
    
## 6/2/14
### Literature Review
* ***Fluidigm SNPtrace Panel***
    * 96 SNP panel used to detect sample contamination and distinguish related samples
    * 6 gender SNPs divided evenly between X and Y chromosomes used to identify contamination
* ***Detecting and Estimating Contamination of Human DNA Samples in Sequencing and Array-Based Genotype Data***
    * Likelihood-based methods with sequence data and array-based genotype data
        * maximizing likelihood of base reading for a value between 0 and 1 of contamination level
        * set of genotypes for each sequence sample is known and authors investigated whether sequencing reads all originated from the targeted sample
    * Likelihood-based methods using just sequence data
        * similar method but prior genotype data unknown
    * Likelihood-based methods just using array-based data
        * used Illumina Infinium florecense
        * Guassian distribution of A and b allele intensity data
        * Maximized likelihood using a grid search on interval [0,1/2] followed by Brent's algorithm
    * Regression-Based method using just array-based data
        * Linear regression of BAF (B allele frequency)
    * Compared their study to heterozyogosity and HET/HOM ratio used to find contaminated DNA in type 2 diabetes study
        * found that the results of methods were consistent but likelihood-based methods and regression-based methods were more sensitive than heterozygosity-based methods
* ***A Novel Method for Detecting Contaminated Samples Based on Illumina Sequencing Data***
    * Also mention HET/HOM but state that populations with recent admixture will skew towards heterozygosity while populations with inbreeding will skew towards homozygosity
    * develop Mappability score
    * compare SNP sites on X and Y chromosome to test for contaminationn
    * use data from unique SNP sites to find contamination
    * I think this paper was translated into English; I had a really hard time understanding what they were tyring to say
* ***ContEst: estimating cross-contamination of humnan samples next-generation sequencing data***
    * General population frequency info and sequencing data in BAM format
    * Bayesian approah to calculate the posterior prob. of the contamination level and determine maximum a posteriori probability (MAP) estsimate of contamination level
        * Indentified homozygous SNPs
        * Bayes law to calculate posterior probability
        * Used uniform prior on contaminatinon fraction
        * Quality of bases typically represented using a Phred-like Q-score (auality = probability of read being incorrect)

### Other
* Also did some work trying to mathematically relate likelihood ratio to the heterozygousity proportion.  No significant break through yet.
* Began reading the MCMC notes (read sections 1 and 2)

## 6/3/14
### Some more literature review
* ***Source-Sink Estimates of Genetic Introgression Show Influence of Hatchery Strays on Wild Chum Salmon Populations in Prince William Sound, Alaska***
    * Used microsatellites and Hardy Weinberg Equilibrium for quality control to identify contaminated samples
    * HWE analysis used similar methods as us (i.e. calculated probability of genotypes in contaminated and noncontaminated  individulas) but used MCMC methods rather than likelihood ratio test to identify contaminated samples
    * used 7 microsatellites to identify contamination
    * Results
        * 1.2% of contemporary scale samples tested positive for contamination using HWE but not microsatellites
        * 55.2% of old samples tested positive for contamination, but the microsatellite method identified much more contamination than the HWE method
        
### Learning about MCMC
* Sections 3 - 5
* had some difficulty understanding Gibbs Sampling
    * this site helped me understand better: http://pareto.uab.es/mcreel/IDEA2014/MCMC/mcmc.pdf

### Editing Code
* added input for a title of the ROC graph in the `ROC` function
* added and edited documentation in all functions except `likelihood1
* added examples to `ROC`, `lratio`, and `hetero`

## 6/4/14
### Learning about MCMC
* Finished reading the notes on MCMC