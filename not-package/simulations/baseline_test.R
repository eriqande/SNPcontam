### Look for contaminated samples in the baseline
populations <- split(x = swfsc_chinook_baseline,f = swfsc_chinook_baseline$Pop )
baseline_test <- lapply(populations, function(x) {single_pop_MCMC(data = x, locstart = 5)})
z_results <- lapply(baseline_test, function(x) {x$analysis$z_pm})
contam_z <- lapply(z_results, function(x) {which(x >= 0.5)})
contam_id <- lapply(1:69, function(x) {if(length(contam_z[[x]]) > 0) {populations[[x]]$ID[contam_z[[x]]]}})
z_pm <- lapply(1:69, function(x) {if(length(contam_z[[x]]) >0) {z_results[[x]][contam_z[[x]]]}})
missing_data <- rowSums(is.na(swfsc_chinook_baseline[unlist(contam_id),-(1:4)]))/2
baseline_results <- data.frame(z_pm = unlist(z_pm), missing_data = missing_data)

### Look for contaminated samples in the CA Department Fish and Game genotypes
p_data <- prepare_base_and_mix_for_mixed_MCMC(B = swfsc_chinook_baseline, B_locstart = 5, B_pops = swfsc_chinook_baseline$RepPop, M = cdfg_port_samples, M_locstart = 1)
p_contam <- Pcontam(snp_zeroes = p_data$zeros, snp_ones = p_data$ones, genos = p_data$mixmat, lambda = .5)
p_clean <- P_likelihood(snp_zeroes = p_data$zeros, snp_ones = p_data$ones, genos = p_data$mixmat, lambda = .5)
cdfg_MCMC <- mixed_MCMC(data = p_data$mixmat, contam_data = p_contam,clean_data = p_clean, inters = 1000)
MCMC_output <- analyzed_mixed_MCMC(MCMC = cdfg_MCMC, burnin = 100, B_pops = levels(swfsc_chinook_baseline$RepPop))
contam_id2 <- which(MCMC_output$z_pm >= .5)
z_pm2 <- MCMC_output$z_pm[contam_id2]
missing_data2 <- rowSums(is.na(cdfg_MCMC[contam_id2,]))/2
cdfg_df <- data.frame(z_pm = z_pm2, missing_data = missing_data2)
cdfg_results <-list(population_id = MCMC_output$pop_id, cdfg_df = cdfg_df)