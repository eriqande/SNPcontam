library(fullsniplings)

#### Setting up the data from the baseline for the function  ####
N = 100
p = .05
N_c = N*p
swfs <- swfsc_chinook_baseline
rownames(swfs) <- swfs$ID
swfs2 <- swfs
swfs2$Pop <- as.integer(factor(as.character(swfs$Pop), levels = unique(as.character(swfs$Pop))))





#### Pull out Contaminated Genotypes and make a mixture sample with them ####
contam <- sample(1:8031,2 * N_c)
contam1 <- contam[1:N_c]
contam2 <- contam[(N_c+1):(2*N_c)]

a <- sample((1:8031)[-contam], N - N_c)

genos <- swfs2[a,-(1:4)]   # uncontaminated fish in the mixture
geno1 <- get_snp_genos(swfs2[contam1,-(1:4)])$mat  # for contaminated fish in mixture
geno2 <- get_snp_genos(swfs2[contam2,-(1:4)])$mat
geno_c <- geno1 + geno2
geno_c[geno_c == 4] = 2
geno_c[geno_c == 3] = 1
geno_c[geno1 ==1 & geno2 == 1] = 1
geno_c[(geno1 == 0 & geno2 == 2)|(geno1 == 2 & geno2 == 0)] = 1
geno_c[is.na(geno1) & is.na(geno2) == FALSE] = geno2[is.na(geno1) & is.na(geno2) == FALSE]
geno_c[is.na(geno2) & is.na(geno1) == FALSE] = geno1[is.na(geno2) & is.na(geno1) == FALSE]

colnames(geno_c) <- paste(colnames(geno1), colnames(geno2), sep="-x-")

pop_id <- swfs2$Pop[a]
contam_id <- rbind(swfs2$Pop[contam1],swfs2$Pop[contam2])

snp_genos <- get_snp_genos(genos)
data <- cbind(snp_genos$mat,geno_c)

# in the end, we drop those from the baseline too
bline <- swfs2[-c(a, contam), ]


get_zeroes_and_ones <- function(x) {
  y <- x[, -(1:4)]
  snp_genos <- get_snp_genos(y)$mat
  snp_indics <- genos_to_indicators(g = snp_genos)
  geno_counts <- count_genos(snp_indics)
  af <- t(alle_freqs(geno_counts, proportion = F))
}

#### Get zero and one allele counts for each population  ####
pop.list <- split(bline, bline$Pop)
alle.counts.list <- lapply(pop.list, get_zeroes_and_ones)
zeros <- do.call(what = cbind, args = lapply(alle.counts.list, function(x) x[,"0"]))
ones <-  do.call(what = cbind, args = lapply(alle.counts.list, function(x) x[,"1"]))


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
#Rep_us <- sapply(1:(N-N_c), function(x) {swfs2$RepUnit[swfs2$Pop==MCMC_pop[x]][1] == swfs2$RepUnit[swfs2$Pop==true_pop[x]][1]})
pop_df <- data.frame(True_Population = true_pop, MCMC_Population = MCMC_pop) #, Same_RepUnit = Rep_us)

# Proportion of Individuals Placed in Correct RepUnit
correct_p = sum(Rep_us)/(N-N_c)

# Assess Mixture Data
count_pop <- as.numeric(table(as.numeric(swfs2$Pop)))
freqs_pop <- count_pop/sum(count_pop)
mixing <- colMeans(test$mixing)
mix_df <- data.frame(True_mix = freqs_pop, MCMC_mix = mixing, Difference = freqs_pop - mixing)


