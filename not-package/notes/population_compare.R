library(fullsniplings)
N = 1000
p = 0.05
N_c = N*p
contamination = TRUE
inters = 1000
data = swfsc_chinook_baseline

#### Get Mixture and Basline ####
b <- make_mixture(data,N,p)


#### Get the Prepared Matrices and Likelihood Matrices ####
MCMC_mats <- prepare_base_and_mix_for_mixed_MCMC(B = b$bline, B_locstart = 4, B_pops = b$bline$Pop, M = b$mixture, M_locstart = 1)
contam_mat <- Pcontam(snp_zeroes = MCMC_mats$zero, snp_ones = MCMC_mats$ones, genos = MCMC_mats$mixmat, lambda = .5)
clean_mat <- P_likelihood(snp_zeroes = MCMC_mats$zero, snp_ones = MCMC_mats$ones, genos = MCMC_mats$mixmat, lambda = .5)

#### Run the MCMC ####
test <- mixed_MCMC(MCMC_mats$mixmat, contam_mat, clean_mat, inters = inters, contamination = TRUE)

#### Create Tables of Population Assignment ####
P = nrow(clean_mat)

# The table pop_means has the individuals from the mixture in the row and all populations in columns
# The table shows the fraction of times that each individual was assigned to each population
# Only creates table for "clean" data. Contaminated samples ommitted
if (contamination){
clean_pops <- lapply(1:(N-N_c), function(x) {test$pops[,x][2*which(test$z[,x] == 0) - 1]})
tmp_pops <- lapply(clean_pops, function(x) {factor(x,level = 1:P)}) 
freq <- lapply(tmp_pops, function(x) {as.numeric(table(x))})
pop_means <- do.call(what = rbind, args = lapply(freq, function(x) x))/inters
colnames(pop_means) <- levels(b$bline$RepPop)
if (p == 0) {rownames(pop_means) <- rownames(b$mixture)}else {rownames(pop_means) <- rownames(b$mixture)[-((N-N_c + 1):N)]}
}else {
  tmp_pops <- lapply(1:N, function(x) {factor(test$pop[,x],level = 1:P)}) 
  freq <- lapply(tmp_pops, function(x) {as.numeric(table(x))})
  pop_means <- do.call(what = rbind, args = lapply(freq, function(x) x))/inters
  colnames(pop_means) <- levels(b$bline$RepPop)
  rownames(pop_means) <- rownames(b$mixture)
}

# Pop_id shows the ID of each individual from the mixtures, it population ID according to the MCMC
# models, and the fraction of times it was assigned to said population
# Contaminated Samples Omitted
tmp_max <- matrix(sapply(1:(N-N_c), function(x) max(pop_means[x,])),ncol=1)
tmp_id <- sapply(1:(N-N_c), function(x) colnames(pop_means)[which(max(pop_means[x,]) == pop_means[x,])])
pop_id <- cbind(tmp_id,tmp_max)
rownames(pop_id) <- rownames(pop_means)

# Table comparing true and MCMC mixing proportions
count_pop <- table(b$bline$Pop)
freqs_pop <- count_pop/sum(count_pop)
mixing <- colMeans(test$mixing)
difference <- as.numeric(freqs_pop) - mixing
mix_df <- data.frame(True_mix = freqs_pop, MCMC_mix = mixing, Difference = difference)
plot(mix_df$True_mix.Freq,mix_df$MCMC_mix,ylim = c(0,.15),xlim = c(0,.18)); abline(a=0,b=1)
