# This script finds the DNA sequences flanking filtered SNPs and downloads from the UCSC Genome DAS Server
# Example url: "http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr9:425000,426000" 
# NOTE: THIS SCRIPT OUTPUTS FASTA FILES WITH HEADERS THAT ARE NON STANDARD.
# REMEMBER TO RENAME YOUR INPUT/OUTPUT FILE BEFORE RUNNING THIS SCRIPT

# Load the relevant package (allows for writing in fasta format)

library(seqinr)

# Sets the working directory

setwd("/home/jlee631/todd_labwork/t1d_snp/wd")

# Import the bed file into R environment

input <- read.table("PDX1_t1d_snp.bed", sep = "\t", as.is = TRUE)

# The next part sets up all the objects/strings required for the loop
# This sets the base url that we will append everything we'll change onto 
# You can see the reference genome in this url is hg19.

base.url = "http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment="

# Establish the objects for writing the sequences and the names of the sequences to

sequences <- data.frame()
names <- data.frame()

# The for loop will go through each individual line of the .bed file and take 15 bases on each side flanking the SNP

for (i in 1:length(input$V2)) {
  # We are appending "chr#", ":" "start position", "," "end position" to base.url for each line to loop through
  # e.g. "http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr9:425000,426000" 
  # This sets up the three variables that will change while looping
  chr <- input$V1[i]
  start.pos <- input$V2[i]
  end.pos <- input$V3[i]
  # Changes the start/end position to get flanking bases
  start.pos <- start.pos - 15
  end.pos <- end.pos + 15
  
  # Creates the full url by putting together all the defined variables from before
  full.url <- paste(base.url, chr, ":", start.pos, ",", end.pos, sep = "")
  
  # Reads the sequence from the UCSC DAS Server and takes only the sequence portion
  full.text <- readLines(full.url)
  sequence.start <- grep("<DNA length", full.text)
  full.text <- full.text[-1:-sequence.start]
  tail.text <- grep("</DNA>", full.text)
  length.text <- length(full.text)
  full.text <- full.text[-tail.text:-length.text]
  
  # Appends the sequences and names to the established objects
  sequences <- c(sequences, full.text)
  names <- c(names, paste(input$V1[i], input$V4[i], sep = "_"))
}

# Outputs the sequences as a fasta format

write.fasta(sequences, names, file.out = "PDX11_SNP_flanking.fasta", open = "w")
