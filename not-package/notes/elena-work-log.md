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

## 6/5/14
* Read about mixture models

## 6/6/14
* Read more about misture models and MCMC

## 6/9/14
### MCMC Model
* wrote down the full conditional probabilities and the Direct Acyclic Graph (DAG) for an MCMC model to find probability of contamination and allele frequencies
    * Equations and graph written down in my notebook
* wrote code to implement MCMC method
    * using function `contam_MCMC` for algorithm
    * function `full_z` used to calculate full conditional distribution of z

### Plans for Tomorrow
* input some of the equation derivation into the `main-body-text.tex` manuscript
* continue commenting and documenting MCMC equations

## 6/10/14
### Manuscript
* began to write down the derivation of the MCMC equations in `main-body-text.tex`

### MCMC code
* wrote some of the documentation for `contam_MCMC` and `full_z`
  * still need to include examples and some description
  
## 6/11/14
### MCMC code
* changed `full_z` to get rid of the `apply()` command and used `log()`, `colSums()`, and `exp()`. Also changed output from list to just `return(p)`
* changed `contam_MCMC` to reflect change in the `full_z` output
* created `simulate_genos` to simulate random contaminated and clean genotypes with a certain proportion of contamination
* created `test_MCMC` to compare MCMC of simualted samples to true results

##6/12/14
### MCMC simulation
* created `loci_table`
    * compares the results of the MCMC for different number of loci
    * outputs a table
* created `missingl_table`
    * compares the results of the MCMC for different proportion of missing loci
    * outputs a table
* created `threshold_table` 
    * to compare the results of the MCMC for different threshold levels
    * outputs a table

### Things to go over with simulation
* sample with replacement from FRHSP or equivalent \theta (+ unif(-.01,0.1) so there is difference between reps) and sample size @ 200
    * Number of loci (20,60,100,200)
    * \rho (.025, .075, .2, .5)
* store \theta in case we want analysis 
* take the last 900
* do 100 reps with different \theta
* record how long it takes: proctime
* look into the package parallel
    * mclapply()
    
## 6/13/14
* started working on the simulation called `MCMC_sims`
* changed `test_MCMC` to have different outputs that correspond better with the desired simulation output

## 6/16/14
### Lab Work
* practiced pipetting with Vanessa
* went through preparing PreAmp and running PCR on for SNP chip

### Simulation
* used `lapply` to run the `test_MCMC` function on all of the
* set up dataframes with the necessary information for making the three different sets of graphs
* set up function `MCMC_plots` that will create the three graphs
* wrote code to create the boxplots for the posterior mean of z

## 6/17/14
### Lab Work
* went through the preparation of the SNP chip and practiced loading old SNP chip wih water
* learned some of the protocol for PCR

### Simulation
* wrote code for all the graphs in `MCMC_zplots`, `MCMC_alleleplot`, and `MCMC_rhoplot`
    * have two different graphs for `MCMC_zplots` because I couldn't decide which one looks the best
* added graph functions to `MCMC_sims` so that the code now produces graphs

## 6/18/14
### Rccp
* read about Rccp
* started working on code to calculate the likelihood of an individual coming from each population in `P_likelihood`

### Talked with Anthony and Carlos about Lab Experiment
* Carlos is going to send Coho data where contaminated individuals were identified using microsatellites
* Anthony is thinking of more details about an experiment investigating effect of concentration of contaminated and clean DNA

## 6/19/14
### Rccp
* finished Rcpp code
* didn't know what the input data would actually look out, so it might be missing a few steps
    * the inputs are
        * genos: a matrix with genotypes where the number of rows equal to loci and number of columns equal to individuals, i.e the ith column is the genotype of the ith individual
        * gc: a matrix with number of rows equal to number of populations and number of columns equal to the number of loci. The matrix is filled with allele frequencies of the 1 allele at a particular locus in a particular population

### Simulation
* changed code in `MCMC_sims` so that we can change the rho values and number of loci used and still use the same code
* added `allele_table` which makes a table of the mean difference between the allele value and the poterior mean for each of the 16 different scenarios
    * I used xtable, but I am not sure how this will actually turn into a latex table

## 6/23/14
* fixed sum of the `MCMC_sims` code to be more efficient and worked with produced figures and graphs

## 6/24/14
* continued to work on the figures and graphs from the `MCMC_sims` data
* continued to fix some of the `MCMC_sims` code to get rid of loops
    * used `do.call()` to combine data frames rather than setting up large complicated data frame in which to store data
    
## 6/25/14
* wrote a vignette to explain the how to run `MCMC_sims` and the plotting functions
    * also wrote information regarding `MCMC_sims` in the README for the package
* tried to run `MCMC_sims` on the computer with many cores, but there was something wrong with the simulation code

## 6/26/14
* fixed the problem with the simulation code so that genotype pool with 0 contaminated samples can be simulated
* removed all the simulation functions from the R package and put in `simulation_functions.R` in `not-package/simulations`
* made Rmarkdown documents `02_makecharts_1.R` and `03_maketables_1.R` contain the code for making charts and tables

## 6/30/14
### Directed Acyclic Graph
* Created a DAG using Inkscape for the one population MCMC method
* downloaded `Macports` and `pdf2svg` in order to be able to import Latex text from Latexit to Inkscape as an svg

### Rcpp Code
* Changed the input matrices for `P_likelihood` code
    * added two separate matrices for zero allele counts and one allele counts

## 7/1/14
### Rcpp
* created function `Pcontam` which creates a matrix of likelihoods for a contaminated sample
    * the table contains the likelihood that each sample is the result of contamination between all combinations of two populations
* fixed problems with `P_likelihood` function
    * fixed loop so that the likelihood value reset for each individual and population
* added documentation to `P_likelihood`

### MCMC code
* added documenttaion to `analyze_MCMC` function

### DAG
* added the DAG_1 to an inkscape folder in `not-package`
* changed the orientation of the DAG so that it is less narrow

### Mixture Model
* formulated full poterior distributions for the MCMC mixture model

## 7/2/14
### DAG
* Began working on `DAG_2`
    * a DAG for the MCMC mixture model

### Mixture Model
* created function `mixed_MCMC`
    * runs the MCMC for the mixture model
    * updates $\rho$ (contamination proportion), $z_i$ (contamination id), $u_i$ (population or population pair id), and $\pi$ (mixture proprotions)

## 7/3/14
### DAG
* Finished `DAG_2`

### Mixture Model
* Finished `mixed_MCMC`
* Created Rscript `sim_mixture` that randomly picks individuals from the baseline and then test the mixture model, `mixed_MCMC`
    * created `README` in `not-package\notes` to explain the code
    * `sim_mixture` does not simulate contaminated samples
    * after running the code, it seems that the model has a low accuracy (~50%) for assigning individuals to the correct RepUnit
    * with 50 individuals, 1000 sweeps, and 95 loci, `mixed_MCMC` took about 10 seconds
    
## 7/7/14
### Mixture Model
* Added contaminated samples to the `sim_mixture` script
* Went over the `sim_mixture` script with Eric and made some changes to simplify code
* Looked at the GSI of data to see if there might be a bug in `mixed_MCMC`
* New functions:
    * `make_mixture` makes a baseline and a mixture matrix for a simulation when given the original baseline
    * `get_likelihood_matrices` computes the likelihood matrices for clean and contam given a baseline and a mixture matrix
* Also made script `population_compare` which runs `make_mixture`, `sim_mixture`, and `mixed_MCMC` and displays population assignment info for non-contaminated individuals
    * one table contains the fraction of times that each individual was assigned to each population
    * the other contains the ID of the individual, the assigned population using the MCMC, and the fraction of times it was assigned to the population

## 7/8/14
### Mixture Model
* Worked on trying to figure out why MCMC is bad at population assignment

#### Test 1: LKuskokwimBristolBay, CentralValleysp--Butte_Cr_Sp, and NPugetSound
* All CentralValleysp fish identified correctly
* All NPugetSound fish identified although all identified to Marblemount and not to separate pops within NPugetSound
* No LKuskokwimBristolBay identified correctly -- all identified to Marblemount
* When run with only CentralValleysp and LKuskokwimBristolBay, only one fish was put into LK
* the likelihood values for LK are consistently lower than NPugetSound and CentralValleysp, which probably contributes to the problem with identification, but I checked the likelihood values and they are correct
* Marblemount also has higher likelihood compared to the other NPugetSound population

#### Test 2: CentralValleysp--Mill_Cr_sp, CentralValleysp--Deer_Cr_sp, NOregonCoast, NSEAlaskaChikatR, SSEAlaska
* All CentralValleysp identified, but all but one to Mill_Cr_sp
* All NOregonCoast identified, but all but one to Alsea
* No Alaskan fish id, all identified as NOregon Coast fish
* running this test with $\rho = 0$ did not change the population assignment
* also tested both Test 1 and Test 2 with the prior for $\pi$ equal to a Dirichlet distribution with 1 as each population's parameter and it did not change results significantly

#### Ran on Whole Baseline without Contamination
* did a similar job of IDing populations, and still gets nothing in Alaska
* the fraction of population assignments to the max population is much lower than with contamination parameter
* did a much better job of estimating $\pi$

#### Found the mistake!!!!!
* I was not updating $\pi$ correctly 
* now the MCMC seems to run better and the values for $\pi$ are much more reasonable
* still can't identify populations in Alaska

#### Code
* Added to the `mixed_MCMC` code so that it can be run with $\rho = 0$ and only run the MCMC for mixing proportions and population assignment
* changed `make_mixture` so makes different mixture matrix and baseline if $\rho = 0$
* changed the `population_compare` script to be general so that it can be run for any portion of the baseline
    * also changed code so that it can be run with $\rho = 0$

### Presentation
* Created rough outline for the presentation using the the LaTex template
* Looked up parameters for the talk: 15 minutes **including** questions (so really only 12 minutes)

## 7/9/14
#### Presentation
* Added pictures and text to the presentation slides
* started working on the mathematics slide, which is unorganized right now

## 7/10/14
#### Presentation
* Added more to the mathematics slide and added a few more sentences to other slides
#### Mixture Model
* changed `make_mixture` so that lists of IDs can be input for both contaminated and clean individuals

## 7/11/14
#### Presentation
* began working on abstract
#### Mixture Model
* changed `make_mixture` so that mixture matrix and baseline matrix have the same format (pre-transformation to 0s,1s, and 2s)
* Eric fixed the mistakes with the code that was counting the 0 and 1 alleles
    * I retested the MCMC, and it is much better at identifying populaitons (especially Alaska)
    * mixing proportions are very accurate
    
## 7/14/14
### Presentation
* finished and sent abstract to the NOAA student scholarship people

### Mixture Model
* changed the `make_mixture` code so that it includes RepPop, RepUnit, and Pop in the mixture and baseline
* changed the flow in `population_compare` so that it works with `prepare_base_and_mixe_for_MCMC`
* also added code in `population_compare` to calculate the fraction of samples that are correctly assigned to their RepUnit and population
* began working on simulation for the mixture model
    * created `mixed_simulation_function` to store all of the simulation functions and moved `make_mixture` into the script
    * created code `make_lists` to take the proportions of the from the `ca_fishery_proportions` and take individuals to be used in the mixture

## 7/15/14
### MCMC Simulation
* finished code for the MCMC simulation in `mixed_MCMC_sims`
    * created function `mixed_MCMC_with_means` to calculate posterior means and such so that an lapply can be applied in `mixed_MCMC_sims` to run all the simulations
    * briefly checked the MCMC and it seems to be functioning correctly
* still need the make the graphs and tables

## 7/16/14
### MCMC Simulation
* added code to the simulation function `mixed_MCMC_sims` so that it uses different rho values
    * now code runs the MCMC for all combinations of fisheries proporitons and rho values
* added code to `mixed_MCMC_functions` to make graphs and tables and to calculate fraction of individuals correctly identified
* created scripts to run the simulation and to make the graphs and tables for the mixed MCMC
* created script to run both tables and make graphs and tables

## 7/17/14
### Presentation
* added images of salmon
* added figures from the mixed_MCMC function
* wrote text for the mixed_MCMC simulation and the conclusion
* now drafts of all slides are finished

## 7/18/14
### Mixture Model
* changed the tables and graphs
    * changed z-table so that the fisheries are horizontal and the rho values are the columns because the table was to big in the other direction
    * took off the x-axis labels and tick marks for the rho plot

### Single Pop MCMC model
* added function `single_pop_MCMC` that wraps all of the MCMC functions for the single population (`analyze_MCMC` and `contam_MCMC`)
    * also included function within `single_pop_MCMC` called `change_to_ones_zeros_twos` that changes the data (which is in the 2 column SNP format) to data with individuals in the columns and loci on the rows with only zeros, ones, twos, and NAs
    * made `change_to_ones_zeros_twos` from code in `prepare_base_and_mix_for_mixture_MCMC` but edited for only one set of data instead of mixture and baseline
* also cleaned up the R directory because there were unused functions from before when I was testing the MCMC model

### Presentation
* practiced presentation with Eric and wrote down some changes to make to my presentation slides and presentation

## 7/21/14
### Presentation
* made some changes to the presentation slides
    * added urls for all the photos that I used
    * took out the allele CI table in the single population results
    * changed the test in the Conclusion slide
    * added GitHub respository information to the last Mathematics slide
* made changes to my presentation text that Eric and I talked about on Friday and also tried a new way to explain Bayesian Statistics (I am going to try it out on my friend later)

### Coho Data
* loaded the Coho data into R and saved three separate files
    * `kl_genotypes` are the samples from the Klamath watershed
    * `noyo_genotypes` are the samples from the Noyo River
    * `contaminated_coho` are the IDs, origins, and contamination status of all the suspect samples
* ran the Klamath data through the single population model: `single_pop_MCMC`
    * seemed that the model was not doing a very good job of finding contaminated individuals
    * there it found 6 of 36 of the Klamath samples that were "cross contaminated" and identified non fo the "possibly contaminated" samples (only one microsatellite locus had more than 2 alleles)
    * also seemed like it was only identifying the contaminated individuals with lots of missing data
* have not tried to test and mixture model type stuff yet

## 7/22/14
### Presentation
* took out a the slide about heterozygousity because it was too long, and I think it is not necessary for the presentation
* also took out the NMFS endangered populations graph because explaining about it was a little to long
* timed presentation and it was about 14 minutes, so I need to shorten it by 2 minutes
* changed some minor typos
* Presented to lab:
    * total presentation was about 17 minutes, but people were asking questions and giving suggestions throughout the presentation
    * suggestions to change minor typos
    * also suggestions to add slide about contamination and more clearly explain the origin of contamination
    * suggestions to add an acknowledgement slide

### Simulations
* added documentation to some of the functions that were missing documentation
* re-added the noyo river data to the not-package/data folder because I had missentered it yesterday
    * there are 7 individuals that are cross-contaminated, possibly contaminated, or miss indentified species?, and the MCMC is not seeming to find any of them

## 7/23/14
### Presentation
* added the suggestions from the lab presentation on Tuesday
    * main changes: added contamination slide, added acknowledgment slide 
* also changed some of my talk

## 7/24/14
### Manuscript
* added some text to the "Methods" section of the manuscript describing the second model
* began to work on writing the derivation down in the appendix for the second model

### Model Simulations
* ran the baseline data through the single population MCMC and ran the CDFG samples through the mixture model MCMC
    * mostly picked up samples with a lot of missing data, but contamination was picked up in a few samples with little or no missing data
    
## 7/25/14
### Manuscript
* worked on the derivation in the appendix for the mixture model

### Presentaiton
* made minor changes on the presentation
    * changed the "P" in the probabilities to upper case
    * added a line in the conclusion about next steps
