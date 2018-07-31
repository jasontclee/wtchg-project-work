# 2018-07-31
# This script is to extract the SNP rsid from the last column and match it with the other file loaded.
# It is custom ID friendly, meaning the custom IDs I generated if the SNPs did not call back an rsid in rAggr are compatible with this script (but not others)

# Set working directory

setwd("/home/jlee631/todd_labwork/t1d_snp/3-1_hg19_MAGIC")

# Load relevant libraries

library(tidyr)
library(plyr)

# Load the data into the R environment

agg_snp <- read.table("agg_snp_update.csv", header = FALSE, sep = ",", as.is = TRUE)
proins <- read.table("agg_snp_proins.csv", header = FALSE, sep = ",", as.is = TRUE)

# Separate the last column

agg_snp_sep <- separate(data = agg_snp, col = V4, into = c("rsid", "info"), sep = ":", extra = "merge")
proins_sep <- separate(data = proins, col = V4, into = c("rsid", "info"), sep = ":", extra = "merge")

overlay <- merge(agg_snp_sep, proins_sep, by.x = "rsid", by.y = "rsid")

#agg_snp_sep_id <- agg_snp_sep$rsid
#proins_sep_id <- proins_sep$rsid

#overlay <- proins_sep_id %in% agg_snp_sep_id

write.table(overlay, file = "overlay.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
