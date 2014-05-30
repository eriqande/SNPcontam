---
title: "Elena's Worklog"
author: "Elena Venable"
output:
  html_document:
    toc: true
    theme: united
  pdf_document:-
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
    + Wrote the documentation for `lratio`, `random_gene` 