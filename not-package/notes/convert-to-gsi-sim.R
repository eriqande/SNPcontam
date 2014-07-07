
# https://github.com/eriqande/gpiper.git
library(gpiper)
library(SNPcontam)

# convert the baseline to gsi sim format
aa <- swfsc_chinook_baseline
rownames(aa) <- aa$ID
aa <- aa[, -(1:4)]
aa[is.na(aa)] <- 0
gPdf2gsi.sim(d = aa, pop.ize.them = factor(swfsc_chinook_baseline$RepPop, levels=unique(swfsc_chinook_baseline$RepPop)))


# to make just a single file for the mixture from a subset of individuals you could do like this
bb <- aa[sample(1:8031, 100), ]
gPdf2gsi.sim(d = bb, outfile = "fakemix.txt")



# could run gsi_sim like this
#   gsi_sim -b gsi_sim_file.txt -t fakemix.txt  


