# This R script will function to take the raw csv files from rAggr and format them for eventual analysis
# We will process the file in .csv format as it is given, remove unnecessary columns (only retain chromosome, position, "name", and r^2 score)
# A new column will also be added for the position to make it consistent with .bed format (the position that flanks the SNP will be "-1" what is stated e.g. pos 2666 with be 2665-2666)
# Finally, we will convert the file to bed format

# First, load the packages

library(tidyr)

# Set the working directory

setwd("/home/jlee631/todd_labwork/t1d_snp/raw T1D SNP data/410012392_union")

# Input the raw csv file into the working environment

agg_snp <- read.csv("410012392.csv", stringsAsFactors = FALSE)


# Have to do a bit of workaround to make this safer for the future
# I'm first going to convert anything in the X chromosome as chr 24, Y chromosome as chr 25
# This will allow us to coerce the chr column as a number, otherwise the sorting gets messed up
# I will then convert it back to a string later on

agg_snp$SNP2.Chr <- sub("X", "24", agg_snp$SNP2.Chr)
agg_snp$SNP2.Chr <- sub("Y", "25", agg_snp$SNP2.Chr)

# This will now force a conversion from strings to numeric. If you get an error here, it means there was a chromosome that was not 1-23, X or Y (so you should check the raw data)

agg_snp$SNP2.Chr <- as.numeric(agg_snp$SNP2.Chr)

# This deletes all the columns we do not want

agg_snp <- subset(agg_snp, select = -c(SNP1.Pos, SNP1.Ref, SNP1.Alt, SNP1.Minor.Allele, SNP2.Ref, SNP2.Alt, D., SNP1.Chr, SNP2.Minor.Allele))

# We are going to merge all the relevant information into a single column
# SNP 2 - SNP2MAF - SNP 1 - SNP1MAF - Distance - Population (more description given in the final documentation at the end)

agg_snp <- unite(agg_snp, SNP2.Name, c(SNP2.Name, SNP2.MAF), remove = TRUE)
agg_snp <- unite(agg_snp, SNP2.Name, c(SNP2.Name, SNP1.Name), remove = TRUE)
agg_snp <- unite(agg_snp, SNP2.Name, c(SNP2.Name, SNP1.MAF), remove = TRUE)
agg_snp <- unite(agg_snp, SNP2.Name, c(SNP2.Name, Distance), remove = TRUE)
agg_snp <- unite(agg_snp, SNP2.Name, c(SNP2.Name, Population), remove = TRUE)

# This will order the SNPs to set up duplicate removal

agg_snp <- agg_snp[order(agg_snp$SNP2.Pos, -agg_snp$R.squared),]

# This will now delete all the duplicates and leave the SNPs with the highest R squared to the index SNPs
# Note: It seems like 2 index SNPs can map to a linked SNP, this deletion step will delete the one that is lower in R squared

agg_snp <- agg_snp[!duplicated(agg_snp[c("SNP2.Pos")]),]

# Now we are going to sort these SNPs back into order by chromosome

agg_snp <- agg_snp[order(agg_snp$SNP2.Chr),]

# Add the starting position for the SNP

agg_snp$start.pos <- agg_snp$SNP2.Pos - 1

# Converting chromosome 24 and 25 back to X and Y

agg_snp$SNP2.Chr <- sub("24", "X", agg_snp$SNP2.Chr)
agg_snp$SNP2.Chr <- sub("25", "Y", agg_snp$SNP2.Chr)

# Need to add the string "chr" to the front of the chromosome numbers

agg_snp$SNP2.Chr <- paste("chr", agg_snp$SNP2.Chr, sep = "")

# This reorders the data frame to resemble the bed format (without the start position column)

agg_snp <- agg_snp[c("SNP2.Chr", "start.pos", "SNP2.Pos","SNP2.Name","R.squared")]

# Nowe we can save the file
# The file columns are as follows:
# Chromosome number (e.g. "chr6")
# SNP starting position (e.g. 1200)
# SNP ending position (e.g. 1201)
# The name of the row, the following information is delimited by "_": 
#          "SNP 2 ID", "SNP 2 MAF", "Index SNP ID", "Index SNP MAF (Minor allele frequency)", "Distance between Index and Linked SNP", "Population from 1000G Project"
# R^2 value to the population listed

write.table(agg_snp, file = "agg_snp.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(agg_snp, file = "agg_snp.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)






