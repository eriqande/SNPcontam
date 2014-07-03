# Testing the Mixture MCMC with Non-Contaminated Data


## Specify the Desired Number of Individuals
In this case, I choose to randomly select 50 individuals from the baseline so:
```
N = 50
```

## Prepare the Data Set
Using the code below, load the baseline data and change the population names to numbers.  In this case, the populations should now be numbers 1 through 69 because this data set contains 69 different populations.
```
swfs3 <- swfsc_chinook_baseline
swfs4 <- data.frame(lapply(swfs3, as.character), stringsAsFactors=FALSE)
Pop <- swfs3$Pop
swfs4$Pop[1] = 1
k = 1; 
for (i in 2:8031){ 
  if (swfs3$Pop[i-1] != swfs3$Pop[i]){
    k = k+1}
  swfs4$Pop[i] = k
}
```

## Select the Individuals from the Baseline
The code below randomly samples the individuals from the baseline, stores their population identification (a number 1 through 69), and trims first two columns of of the data so that it can be passed the functions for receiving genotype information.

```
a <- sample(1:8031,N)
pop_id <- swfs4$Pop[a]
genos <- swfs4[a,-(1:2)]
```

## Change the Genotype Data into Suitable Data for the MCMC Code
The function`get_snp_genos` creates a matrix with number of columns equal to the number individuals (N) and number of rows equal to the number of loci.  The genotypes of the individuals are represented as counts of the 1 allele and are values of 0, 1, or 2.

```
snp_genos <- get_snp_genos(genos)
data <- snp_genos$mat
```

## Get the Allele Counts for Each Population at Each Locus

The code below is used to find the allele counts, which are needed for the functions `P_likelihood` and `Pcontam`.

```
# Creates a list where the data from each population is a separate list item
boing <- lapply(1:69, function(x) {subset(swfs4,swfs4$Pop == x)[,-(1:2)]})

# Gets the genotype counts for each population
tmp <- lapply(boing, function(x) {snp_genos = get_snp_genos(x)})
snp_indics <- lapply(1:69, function(x) {genos_to_indicators(g = tmp[[x]]$mat)})
geno_counts <- lapply(snp_indics, function(x) {count_genos(x)})

# Creates the zero matrix, which is the count of zero alleles at each loci for each population
zeros_t <- lapply(geno_counts, function(x) {x[1,]*2 + x[2,]})
zeros <- matrix(unlist(zeros_t),ncol=69)
# Creates the one matrix, which is the count of one alleles at each loci for each population
ones_t <- lapply(geno_counts, function(x) {x[3,]*2 + x[2,]})
ones <- matrix(unlist(ones_t), ncol = 69)
```

## Calculate the Likelihood Matrices
The code below is used to calculate the likelihood matrices for contaminated and non-contaminated assumptions.
```
clean_prob <- P_likelihood(zeros,ones,data,.5)
contam_prob <- Pcontam(zeros,ones,data,.5)
```

## Run the MCMC
The MCMC data is calculated using `mixed_MCMC`.
```
test <- mixed_MCMC(data, contam_prob, clean_prob, inters = 1000) # took about 9 seconds with N = 50
```

## Calculate and Compare Stats
Below are some options for computing and displaying some of the data produced by `mixed_MCMC`.  

```
# Posterior Mean of the Contam Proportion (rho)
rho <- mean(test$prob_contam)

# Table for Population Data
pops <- 1:69
true_pop <- pop_id # true population id
clean_pops <- lapply(1:N, function(x) {test$pops[,x][2*which(test$z[,x] == 0) - 1]})
tmp_pops <- lapply(clean_pops, function(x) {factor(x,level = pops)}) 
freq <- lapply(tmp_pops, function(x) {as.numeric(table(x))})
MCMC_pop <- sapply(freq, function(x) {pops[which(x == max(x))]}) # MCMC identified pop
Rep_us <- sapply(1:N, function(x) {swfs4$RepUnit[swfs4$Pop==MCMC_pop[x]][1] == swfs4$RepUnit[swfs4$Pop==true_pop[x]][1]}) # logical vector recording if individual was assigned to correct RepUnit
pop_df <- data.frame(True_Population = true_pop, MCMC_Population = MCMC_pop, Same_RepUnit = Rep_us)

# Proportion of Individuals Placed in Correct RepUnit
correct_p = sum(Rep_us)/N

# Table Mixture Data
count_pop <- as.numeric(table(as.numeric(swfs4$Pop)))
freqs_pop <- count_pop/sum(count_pop) # actual mixture frequencies from baseline
mixing <- colMeans(test$mixing) # posterior means mixture frequencies from the MCMC
mix_df <- data.frame(True_mix = freqs_pop, MCMC_mix = mixing, Difference = freqs_pop - mixing)
```