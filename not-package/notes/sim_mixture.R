library(fullsniplings)

# Setting up the data from the baseline for the function
N = 50
swfs4 <- data.frame(lapply(swfs3, as.character), stringsAsFactors=FALSE)
Pop <- swfs3$Pop
swfs4$Pop[1] = 1
k = 1; 
for (i in 2:8031){ 
  if (swfs3$Pop[i-1] != swfs3$Pop[i]){
    k = k+1}
  swfs4$Pop[i] = k
}
a <- sample(1:8031,N)
genos <- swfs4[a,-(1:2)]
pop_id <- swfs4$Pop[a]

snp_genos <- get_snp_genos(genos)
data <- snp_genos$mat

boing <- lapply(1:69, function(x) {subset(swfs4,swfs4$Pop == x)[,-(1:2)]})
tmp <- lapply(boing, function(x) {snp_genos = get_snp_genos(x)})
snp_indics <- lapply(1:69, function(x) {genos_to_indicators(g = tmp[[x]]$mat)})
geno_counts <- lapply(snp_indics, function(x) {count_genos(x)})
zeros_t <- lapply(geno_counts, function(x) {x[1,]*2 + x[2,]})
zeros <- matrix(unlist(zeros_t),ncol=69)
ones_t <- lapply(geno_counts, function(x) {x[3,]*2 + x[2,]})
ones <- matrix(unlist(ones_t), ncol = 69)

# Getting the two probability matrices
clean_prob <- P_likelihood(zeros,ones,data,.5)
contam_prob <- Pcontam(zeros,ones,data,.5)

# Running the MCMC
test <- mixed_MCMC(data, contam_prob, clean_prob, inters = 1000) # took about 9 seconds with N = 50

pops <- 1:69
# Mean contam proportion
rho <- mean(test$prob_contam)

# Create Table for Population Data
true_pop <- pop_id
clean_pops <- lapply(1:N, function(x) {test$pops[,x][2*which(test$z[,x] == 0) - 1]})
tmp_pops <- lapply(clean_pops, function(x) {factor(x,level = pops)}) 
freq <- lapply(tmp_pops, function(x) {as.numeric(table(x))})
MCMC_pop <- sapply(freq, function(x) {pops[which(x == max(x))]}) 
Rep_us <- sapply(1:N, function(x) {swfs4$RepUnit[swfs4$Pop==MCMC_pop[x]][1] == swfs4$RepUnit[swfs4$Pop==true_pop[x]][1]})
pop_df <- data.frame(True_Population = true_pop, MCMC_Population = MCMC_pop, Same_RepUnit = Rep_us)

# Proportion of Individuals Placed in Correct RepUnit
correct_p = sum(Rep_us)/N

# Assess Mixture Data
count_pop <- as.numeric(table(as.numeric(swfs4$Pop)))
freqs_pop <- count_pop/sum(count_pop)
mixing <- colMeans(test$mixing)
mix_df <- data.frame(True_mix = freqs_pop, MCMC_mix = mixing, Difference = freqs_pop - mixing)

