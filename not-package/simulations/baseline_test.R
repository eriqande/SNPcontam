populations <- split(x = swfsc_chinook_baseline,f = swfsc_chinook_baseline$Pop )
baseline_test <- lapply(populations, function(x) {single_pop_MCMC(data = x, locstart = 5)})
z_results <- lapply(baseline_test, function(x) {x$analysis$z_pm})
contam_z <- lapply(z_results, function(x) {which(x >= 0.5)})
contam_id <- lapply(1:69, function(x) {if(length(contam_z[[x]]) > 0) {populations[[x]]$ID[contam_z[[x]]]}})
z_pm <- lapply(1:69, function(x) {if(length(contam_z[[x]]) >0) {z_results[[x]][contam_z[[x]]]}})
missing_data <- rowSums(is.na(swfsc_chinook_baseline[unlist(contam_id),-(1:4)]))/2
baseline_results <- data.frame(z_pm = unlist(z_pm), missing_data = missing_data)
