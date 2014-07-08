
#### Get Mixture and Basline ####
b <- make_mixture(swfsc_chinook_baseline,100,.05)


#### Get the Likelihood Matrices ####
likelihood <- get_likelihood_matrices(b$bline,b$mixture)

#### Run the MCMC ####
test <- mixed_MCMC(b$mixture, likelihood$contam_prob, likelihood$clean_prob, inters = 1000)

#### Create Tables of Population Assignment ####
N = 100
p = 0.5
N_c = N*p

# The table pop_means has the individuals from the mixture in the row and all populations in columns
# The table shows the fraction of times that each individual was assigned to each population
# Only creates table for "clean" data. Contaminated samples ommitted
clean_pops <- lapply(1:(N-N_c), function(x) {test$pops[,x][2*which(test$z[,x] == 0) - 1]})
tmp_pops <- lapply(clean_pops, function(x) {factor(x,level = 1:69)}) 
freq <- lapply(tmp_pops, function(x) {as.numeric(table(x))})
pop_means <- do.call(what = rbind, args = lapply(freq, function(x) x))/1000
colnames(pop_means) <- levels(b$bline$RepPop)
rownames(pop_means) <- colnames(b$mixture)[-((N-N_c + 1):N)]

# Pop_id shows the ID of each individual from the mixtures, it population ID according to the MCMC
# models, and the fraction of times it was assigned to said population
# Contaminated Samples Omitted
tmp_max <- matrix(sapply(1:(N-N_c), function(x) max(pop_means[x,])),ncol=1)
tmp_id <- sapply(1:(N-N_c), function(x) colnames(pop_means)[which(max(pop_means[x,]) == pop_means[x,])])
pop_id <- cbind(tmp_id,tmp_max)
rownames(pop_id) <- rownames(pop_means)

