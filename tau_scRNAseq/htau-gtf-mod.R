# 2018-08-23

# This is a very short script to amend the entries for human Tau in the human GTF before appending to the mouse GTF
# The chromosome will be set as an arbitrary "30" and the locations will all have to be changed to make sure the first base is at "1"

# Set working directory

setwd("/home/jlee631/ref_genome")

# Load the GTF file for modification

htau_gtf <- read.table("ENSG00000186868.gtf", header = FALSE, sep = "\t", quote="")

# Set the chromosome from 17 to 30

htau_gtf$V1 <- 30

# The position that is originally listed is 45894382 (found using htau_gtf)
# Will need to make that value 1 and then everything else shifted relative to that value(basically, subtract 45894381 from everything)

htau_gtf$V4 <- htau_gtf$V4 - 45894381
htau_gtf$V5 <- htau_gtf$V5 - 45894381

# You can sanity check if you want

htau_gtf$V1
htau_gtf$V4
htau_gtf$V5

# Just to make sure this matches with the number originally entered in the ref genome

max(htau_gtf$V5)

# The number winds up being 133953 (which does match)

# Export the file now

write.table(htau_gtf, file = "ENSG00000186868_mod.gtf", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)