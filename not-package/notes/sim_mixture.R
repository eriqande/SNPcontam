library(fullsniplings)

# Setting up the data from the baseline for the function
N = 100
p = .05
N_c = N*p
swfs <- swfsc_chinook_baseline
swfs2 <- data.frame(lapply(swfs, as.character), stringsAsFactors=FALSE)
Pop <- swfs$Pop
swfs2$Pop[1] = 1
k = 1; 
for (i in 2:8031){ 
  if (swfs$Pop[i-1] != swfs$Pop[i]){
    k = k+1}
  swfs2$Pop[i] = k
}
a <- sample(1:8031,(N - N_c))
genos <- swfs2[a,-(1:4)]

# Contaminated Genotypes
contam1 <- sample(1:8031,N_c)
contam2 <- sample(1:8031,N_c)
geno1 <- get_snp_genos(swfs2[contam1,-(1:4)])$mat
geno2 <- get_snp_genos(swfs2[contam2,-(1:4)])$mat
geno_c <- geno1 + geno2
geno_c[geno_c == 4] = 2
geno_c[geno_c == 3] = 1
geno_c[geno1 ==1 & geno2 == 1] = 1
geno_c[(geno1 == 0 & geno2 == 2)|(geno1 == 2 & geno2 == 0)] = 1
geno_c[is.na(geno1) & is.na(geno2) == FALSE] = geno2[is.na(geno1) & is.na(geno2) == FALSE]
geno_c[is.na(geno2) & is.na(geno1) == FALSE] = geno1[is.na(geno2) & is.na(geno1) == FALSE]

pop_id <- swfs2$Pop[a]
contam_id <- rbind(swfs2$Pop[contam1],swfs2$Pop[contam2])

snp_genos <- get_snp_genos(genos)
data <- cbind(snp_genos$mat,geno_c)

boing <- lapply(1:69, function(x) {subset(swfs2,swfs2$Pop == x)[,-(1:4)]})
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

# Posterior Mean of z values
z_pm <- colMeans(test$z)

# Create Table for Population Data
true_pop <- pop_id
clean_pops <- lapply(1:(N-N_c), function(x) {test$pops[,x][2*which(test$z[,x] == 0) - 1]})
tmp_pops <- lapply(clean_pops, function(x) {factor(x,level = pops)}) 
freq <- lapply(tmp_pops, function(x) {as.numeric(table(x))})
MCMC_pop <- sapply(freq, function(x) {pops[which(x == max(x))]}) 
Rep_us <- sapply(1:(N-N_c), function(x) {swfs2$RepUnit[swfs2$Pop==MCMC_pop[x]][1] == swfs2$RepUnit[swfs2$Pop==true_pop[x]][1]})
pop_df <- data.frame(True_Population = true_pop, MCMC_Population = MCMC_pop, Same_RepUnit = Rep_us)

# Proportion of Individuals Placed in Correct RepUnit
correct_p = sum(Rep_us)/(N-N_c)

# Assess Mixture Data
count_pop <- as.numeric(table(as.numeric(swfs2$Pop)))
freqs_pop <- count_pop/sum(count_pop)
mixing <- colMeans(test$mixing)
mix_df <- data.frame(True_mix = freqs_pop, MCMC_mix = mixing, Difference = freqs_pop - mixing)


