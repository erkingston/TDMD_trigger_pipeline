# This script takes as input the list of loop-length selected sites (2_TDMD_trigger_narrower_looplength.R), 
# And the duplex energies predicted by 3_trigger_miRNA_duplex_energy.py
# And it outputs the top 100 candidate sites for manual annotation with BLAST

##############################################################################
library(tidyr)
library(dplyr)
##############################################################################
# Import the relevant data.frames, and combine them into one merged data.frame
energy.df <- read.table("File output from 3_trigger_miRNA_duplex_energy.py", sep = ',', header = FALSE)
select.df <- read.table("File output from 2_TDMD_trigger_narrower_looplength.R", sep = ',', header = TRUE)
colnames(energy.df) <- c('gene', 'energy', 'structure', 'sequence')
merged.df <- merge(select.df, energy.df[c(2:4)], by = 'sequence')

##############################################################################
# Order the data.frame based on energy so the 'best' sites are at the top of the data.frame
# And write out this ordered data frame for use in the final steps of 5_trigger_conservation_analyses.R
ordered.merged.df <- merged.df[order(merged.df$energy),]
write.csv(ordered.merged.df, "Summary duplex energy prediction data for use with 5_trigger_conservation_analyses.R", quote = FALSE, row.names = FALSE)

##############################################################################
# Print out the top 100 sites for the user-driven conservation analysis step
# The sequences of these top 100 sites will be output in 4 .txt files, each composed of 25 sequences in fasta format
# Each of these .txt files can then be fed into the batch BLAST function on the UCSC genome browser to obtain information about location in the genome
# This genomic location information will then be used in XXX to pull out information about the conservation of this region

# First, generate a new data.frame where the sequence data is formatted in fasta file format
reorg.df <- data.frame(information = character())
for (i in 1:100)
{
  name <- paste('>', i, sep = ' ')
  seq <- as.character(ordered.merged.df$sequence[i])
  temp.df <- data.frame(information = c(name, seq))
  reorg.df <- rbind(reorg.df, temp.df)
}

# Then write this fasta file formatted data out to 4 .txt files
for (iter in 1:4)
{
  filename <- paste('Out file path', iter, '.txt', sep = '')
  write.table(reorg.df[(50*(iter-1) + 1):(50*iter),], filename, row.names = FALSE, col.names = FALSE, quote = FALSE)
}