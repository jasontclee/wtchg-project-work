# 2018-07-31
# This script is to extract relevant data from the MAGIC dataset (SNPs that affect proinsulin expression)

# Set working directory

setwd("/home/jlee631/todd_labwork/t1d_snp/3-1_hg19_MAGIC")

# Load relevant libraries

library(plyr)

# Load the data into the R environment

magic_proins <- read.table("MAGIC_proinsulin_for_release_HMrel27.txt", header = TRUE, sep = "\t", as.is = T)

# Data is filtered for a p value < 5e-8 and for positive effect (increases proinsulin expression)

magic_proins_filt <- subset(magic_proins, magic_proins$pvalue < 0.0000001)
magic_proins_filt <- subset(magic_proins_filt, magic_proins_filt$effect > 0)

# Write the data as a list of SNPs

write.table(magic_proins_filt$snp, file = "proinsulin_SNPs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)