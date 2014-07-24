# This is a script for running the mixture MCMC model on the CA Department of Fish and Game
# samples in California ports.
# The output is a list that includes two data frames.  "population_id" is a data frame
# that displays the population assignment of all individuals.  "cdfg_df" is a data frame
# that displays the index of contaminated individuals, their z-value, the number of missing
# loci in the sample, the first pair of populations to which they were assigned, and the second pair of 
# populations to which they were assigned.

library(gtools)
# Prepare the mixture for the MCMC and calculate the one and zero allele counts for each population
p_data <- prepare_base_and_mix_for_mixed_MCMC(B = swfsc_chinook_baseline, B_locstart = 5, B_pops = swfsc_chinook_baseline$RepPop, M = cdfg_port_samples, M_locstart = 1)

# Calculate the likelihood matrices for both clean and contaminated assumptions
p_contam <- Pcontam(snp_zeroes = p_data$zeros, snp_ones = p_data$ones, genos = p_data$mixmat, lambda = .5)
p_clean <- P_likelihood(snp_zeroes = p_data$zeros, snp_ones = p_data$ones, genos = p_data$mixmat, lambda = .5)

# Run the MCMC on the cdfg data
cdfg_MCMC <- mixed_MCMC(data = p_data$mixmat, contam_data = p_contam,clean_data = p_clean, inters = 1000)

# Analyze output and calculate poterior means
MCMC_output <- analyzed_mixed_MCMC(MCMC = cdfg_MCMC, burnin = 100, B_pops = levels(swfsc_chinook_baseline$RepPop))

# Id contaminated individuals and find their z values
contam_id2 <- which(MCMC_output$z_pm >= .5)
z_pm2 <- MCMC_output$z_pm[contam_id2]
missing_data2 <- rowSums(is.na(cdfg_port_samples[contam_id2,]))/2
contam_pops <- MCMC_output$contam_df[,-1]
cdfg_df <- data.frame(z_pm = z_pm2, missing_data = missing_data2, pairs = contam_pops)

# List Output
cdfg_results <-list(population_id = MCMC_output$pop_id, cdfg_df = cdfg_df)
