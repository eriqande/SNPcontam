# This is a script for running the single population MCMC function on the baseline
# data to test for contamination.  
# The output is a data frame including the posterior mean of the z_values of individuals
# identified as contaminated, the ID of the individual, and the amount of missing data
# in the sample

# Split populations so that each population can be tested separatesly
populations <- split(x = swfsc_chinook_baseline,f = swfsc_chinook_baseline$Pop )

# Run MCMC on all 69 populaitons
baseline_test <- lapply(populations, function(x) {single_pop_MCMC(data = x, locstart = 5)})

# Find which samples had z values above 0.5
z_results <- lapply(baseline_test, function(x) {x$analysis$z_pm})
contam_z <- lapply(z_results, function(x) {which(x >= 0.5)})
contam_id <- lapply(1:69, function(x) {if(length(contam_z[[x]]) > 0) {populations[[x]]$ID[contam_z[[x]]]}})
z_pm <- lapply(1:69, function(x) {if(length(contam_z[[x]]) >0) {z_results[[x]][contam_z[[x]]]}})

# Calculate number of missing loci
missing_data <- rowSums(is.na(swfsc_chinook_baseline[unlist(contam_id),-(1:4)]))/2

# Output data frame
baseline_results <- data.frame(z_pm = unlist(z_pm), missing_data = missing_data)
