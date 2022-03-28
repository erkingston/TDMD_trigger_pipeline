# This script parses the output from the user-driven BLAST analysis, 
# And then determines conservation information for all fo the top 100 most energetically favorable sites
# That were output by 4_TDMD_trigger_narrower_energy.R
# The output is the final file that is then used to identify putative trigger sites for experimental validation

##############################################################################
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)
library(plotly)
#############################################################################
# First merge all the various blast information with the sequence information
# Note that as this is written, this is not batchable, and thus each BLAST file must be combined with the corresponding sequence file
# And then the code is manually changed for the next pair and run through again
blast.list <- read.table("Path for BLAST output for sequence set X", sep = '', header = FALSE)
target.list <- read.table("Path for sequence info input to BLAST for set X (1 of the 4 files generated from 4_TDMD_trigger_narrower_energy.R)", sep = ',', header = FALSE)

# Rearrange the target.list, pulling out first a vector of the indices associated with each sequence, and then a vector of the sequences
# Also remove the '>' symbol from the sequence index number
target.df <- data.frame(index = target.list[c(TRUE,FALSE),], sequence = target.list[c(FALSE,TRUE),])
target.df$seq.number <- as.integer(unlist(lapply(lapply(target.df$index, function(x) strsplit(as.character(x), " ")), function(x) x[[1]][2])))

# Select only those regions that shared 100% identity with the sequence queried in the BLAST search
# These regions are identified by the length of the BLAST match (blast.list$V4 == blast.list$V7, and by the 100.0% identity in blast.list$V8)
# Note that this code is required because BLAST will spit out regions of lower identity, or of shorter identity, as well
# Note also that some sequences will blast perfectly to multiple places within the genome, and thus these sequences will have multiple entries
select.blast <- blast.list[(blast.list$V4 == blast.list$V7) & (blast.list$V8 == '100.0%'),c(3,9:12)]

# Occassionally, some of the sequences will fail to blast with 100% identity for the entire queried region
# When that is the case, instead take the blasted region that has the best output score (column V4)
# First check to see whether or not any of the sequence indices are missing from this final select.blast file
# Then for those missing ids, pull out the top scored blasted region
missing.ids <- seq(25)[!(seq(25) %in% unique(select.blast$V3))]
for (id in missing.ids)
{
	sub.df <- blast.list[blast.list$V3 == id,]
	
	if (nrow(sub.df) > 0)
	{
		ordered.df <- sub.df[order(-sub.df$V4),c(3,9:12)]
		select.blast <- rbind(select.blast, ordered.df[1,])
	}
}

# Finally, merge the select.blast and target.df data.frames by the sequence index
colnames(select.blast) <- c('seq.number','chromosome', 'strand', 'start', 'end')
sequence.blast <- merge(target.df[,2:3], select.blast, by = 'seq.number')
write.csv(sequence.blast, "Output file path", quote = FALSE, row.names = FALSE)

#############################################################################
# This code simply takes the files generated from the previous section (which should be 4 files, one for each set of sequences BLASTED on the UCSC genome browser)
blast1 <- read.table("Output file path 1 from the above section", sep = ',', header = TRUE)
blast2 <- read.table("Output file path 2 from the above section", sep = ',', header = TRUE)
blast3 <- read.table("Output file path 3 from the above section", sep = ',', header = TRUE)
blast4 <- read.table("Output file path 4 from the above section", sep = ',', header = TRUE)

all.blast <- rbind(blast1[,2:6], blast2[,2:6], blast3[,2:6], blast4[,2:6])
write.csv(all.blast, "Output filepath for the blast summary information", quote = FALSE, row.names = FALSE)

##############################################################################
# Finally the data has been wrangled in such a way that one can pull out conservation information about each of the candidate trigger sites
# To determine conservation of these sites, we will use the positional information gleaned from blast with the phyloP information pulled from flybase
# I use the dm6.phyloP124way.bed file downloaded from flybase, which has nucleotide level conservation values
cons.info <- read.table("PhyloP file path", sep = '\t', header = FALSE) 
sequences <- read.table("File path for summary blast information (generated in the section above)", ',', header = TRUE)
colnames(cons.info) <- c('chromosome', 'start', 'end', 'phastCons')

# This code determines the mean PhyloP score over the entire sequence pulled out as the trigger site
# The conservation analysis part of this pipeline is admittedly not perfect, and would be something to improve in later iterations
mean.PhyloP <- c()
for (entry in 1:nrow(sequences))
{
	# First pull out information about the genomic location of the site
  chrom <- as.character(sequences$chromosome[entry])
	site.start <- sequences$start[entry]
	site.end <- sequences$end[entry]
	
	# Then isolate the conservation information for the relevant chromosome, and then pull out the information for the specific chromosomal location
	select.cons.info <- cons.info[as.character(cons.info$chromosome) == chrom,]
	start.index <- grep(site.start, select.cons.info$start)
	end.index <- grep(site.end, select.cons.info$end)
	
	# If there's no conservation information for either the start or end position, an NA will be output
	# If both start and end indices have conservation information, then determine the mean conservation across this entire trigger region
	# Note that because na.rm is not specified for mean(), this means that any sites that have NA values for the PhyloP scores in between the start and end indices will yield NA mean values for the conservation
	if ((length(start.index) == 0) | (length(end.index) == 0))
	{
		mean.PhyloP <- c(mean.PhyloP, NA)
	}
	else
	{
		mean.PhyloP <- c(mean.PhyloP, mean(select.cons.info$phastCons[start.index:end.index]))
	}
}

# Generate a summary dataframe with the sequences and the mean PhyloP data
conserved.sequences <- data.frame(sequence = sequences$sequence, conservation = mean.PhyloP)
write.csv(merge(sequences, conserved.sequences, by = 'sequence'), "Output file path for summary dataframe", quote = FALSE, row.names = FALSE)

##############################################################################
# Finally, combine the conservation information and the energy information to obtain a final file that can then be used to select trigger candidates for future experiments
cons.data <- read.table("File path for the file output from the above section", sep = ',', header = TRUE)

# Some of the trigger sequences blast perfectly to multiple locations in the genome (as they come from repetitive sequences)
# This leads to multiple conservation values for these regions
# For these sequences, this code determines the region with the highest conservation value and outputs this as the conservation score for the sequence
max.cons <- c()
unique.seqs <- unique(cons.data$sequence)
chromosome <- c()
start <- c()
end <- c()

# Iterate through each sequence, identify the maximum conservation score associated with that sequence, and output the information corresponding to that site
for (seq in unique.seqs)
{
	sub.df <- cons.data[cons.data$sequence == seq,]
	ordered.sub.df <- sub.df[order(-sub.df$conservation),]
	max.cons <- c(max.cons, sub.df$conservation[1])
	chromosome <- c(chromosome, sub.df$chromosome[1])
	start <- c(start, sub.df$start[1])
	end <- c(end, sub.df$end[1])
	
}
cons.df <- data.frame(sequence = unique.seqs, conservation = max.cons, Chrom = chromosome, start = start, end = end)

# Finally, combine this finalized conservation data with the duplex energy prediction data from 4_TDMD_trigger_narrower_energy.R
# And write it out to a final file that can be used to identify candidate sites for further experimental analyses
energy.data <- read.table("Summary duplex energy prediction data for use with 4_TDMD_trigger_narrower_energy.R", sep = ',', header = TRUE)
energy.cons.df <- merge(energy.data, cons.df, by = 'sequence')
write.csv(energy.cons.df, "Final output file path", row.names = FALSE, quote = FALSE)
