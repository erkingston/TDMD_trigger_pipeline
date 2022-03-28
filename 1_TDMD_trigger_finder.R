library(dplyr)
library(tidyr)
library(seqinr)
library(devtools)
library(stringi)
library(Biostrings)

##########################################################################################
# # This code searches a user-provided transcriptome for potential TDMD triggers for a user-provided list of miNRAs #

##########################################################################################
# # First read in all the necessary files # #
# File 1 = List of miR that are subject to TDMD, column 1 = name, column 2 = sequence
miRNA.list <- read.table("Insert file path", sep = "\t", header=FALSE)
colnames(miRNA.list) <- c('names', 'sequence')

# File 2 = RNAseq counts data for cell-type of interest, column 1 = gene, column 2 = counts
RNA.seq <- read.table("Insert file path",sep = "\t", header=FALSE) 
colnames(RNA.seq) <- c('gene', 'counts')
select.RNA.seq <- RNA.seq[RNA.seq$counts > 10,] # Ask for all genes with more than 10 counts

# File 3 = fasta file with sequences for transcript regions of interest
# For flies, I used 'dmel-all-ncRNA-r6.08.fasta.gz' for ncRNA and 'dmel-all-three_prime_UTR-r6.08.fasta' for 3'UTR (downloaded from flybase)
transcripts <- read.fasta("Insert file path", seqtype = 'DNA', as.string = TRUE)

# Then get the annotations from the fasta file
# Note: parameters in this code must be changed for different fasta files because the annotations are structured differently
# I have noted the parameters used for the drosophila fasta files mentioned above
transcript.annots <- unlist(lapply(transcripts, function(x) attr(x, "Annot", exact = FALSE)))
transcript.names <- unlist(lapply(transcripts, function(x) attr(x, "name", exact = FALSE)))
gene.names <- unlist(lapply(strsplit(sub("* parent=","",transcript.annots), split = ';'), function(x) x[10])) # Toggle to x[10] for ncRNA or x[5] for UTR #
transcript.df <- data.frame(gene.ID = gene.names, sequence = unlist(lapply(transcripts, function(x) x[1])), transcript.ID = transcript.names)
rownames(transcript.df) <- c()

# Choose transcripts to look at based on expression, then generate a DNAStringSet out of them
# Note that this code must be modified based on theh file used
# Once again, I have noted what the code should read for the two fasta files I used
selected.transcripts <- unique(transcript.df[as.character(transcript.df$gene.ID) %in% select.RNA.seq$gene,]) # Toggle this to gene.ID or transcript.ID based on ncRNA or UTR

# This next chunk of code loops through the miRNAs provided in the TDMD-sensitive list, and pulls out potential triggers #
miRNA.targets <- list()
for (miRNA in 1:nrow(miRNA.list))
{
	# First, get the supplemental and seed region of the miRNA sequence #
  # Supplemental region defined as nucleotides 13+
  # Seed defined as nucleotides 2-8
	miR.seq <- as.character(miRNA.list$sequence[miRNA])
	supp <- substr(miR.seq, 13, nchar(miR.seq))
	seed <- paste('T', substr(miR.seq, 2, 8), sep = '')
	
	# Take the reverse complement of these sequences, as these will be what are searched for in the transcriptome
	supp.revcomp <- reverseComplement(DNAString(supp))
	seed.revcomp <- reverseComplement(DNAString(seed))
	
	# First generate a list of sequences that have complementarity to the 3' region #
	# As currently written, the script allows up to 2 mismatches to the supplementary region (max.mismatch = 2)
	# Note that G-U wobbles count as a mismatch, as do indels
	# match.df is a data frame with the gene.id, the full sequence of the transcript that has the match, the nucleotide at which this complementarity starts, and the nucleotide at which this complementarity ends
  match.list <- lapply(selected.transcripts$sequence, function(x) matchPattern(supp.revcomp, DNAString(x), max.mismatch = 2, min.mismatch = 0, with.indels = TRUE))
	match.df <- data.frame(gene.ID = rep.int(selected.transcripts$gene.ID, times = lapply(match.list, function(x) length(x))), full.seq = rep.int(selected.transcripts$sequence, times = lapply(match.list, function(x) length(x))), start = unlist(lapply(match.list, function(x) start(x))), end = unlist(lapply(match.list, function(x) end(x))))
	
	# For each of the regions that had a match to the 3' end of the miRNA, look in the 30nt upstream of the 3' end for seed complementarity for any match to the seed region #
	# The 'query region' is this region in which the seed is search for
	# The full region includes the full string of the region complementary to the 3' end plus the 30 nt upstream, will be printed out in the final data frame
	# Start.end.diffs describes the length of the match to the 3' region of the miRNA, which will be used in downstream analyses looking at the loop length of predicted sites
	query.region <- apply(match.df, 1, function(x) DNAString(substr(as.character(x[2]), as.numeric(x[4]), (as.numeric(x[4]) + 30))))
	full.region <- apply(match.df, 1, function(x) substr(as.character(x[2]), as.numeric(x[3]), (as.numeric(x[4]) + 30)))
	start.end.diffs <- as.numeric(match.df[,4]) - as.numeric(match.df[,3])
	
	# Actually search the query region for a seed match, and then make a data frame for those regions that have both 3' and seed complementarity
	# Note, as written up to 2 mismatches are allowed for the seed region (including G-U wobbles)
	# Also note that one region of complementarity to the 3' end can have multiple regions of complementarity to the seed in different frames, and this code spits out all of them
	# Final dataframe has gene.ID, the sequence of the complementary region (seq), the nucleotide at which complementarity to the 3' region ends (threep.end), the nucleotide at which complementarity to the seed starts (seed.start), and the nucleotide at which the seed complementarity ends (seed.end)
	seed.match <- lapply(query.region, function(x) matchPattern(seed.revcomp, x, max.mismatch = 2, min.mismatch = 0))
	full.match.df <- data.frame(gene.ID = rep.int(match.df$gene.ID, times = lapply(seed.match, function(x) length(x))), seq = rep.int(full.region, times = lapply(seed.match, function(x) length(x))), threep.end = rep.int(start.end.diffs, times = lapply(seed.match, function(x) length(x))), seed.start = unlist(lapply(seed.match, function(x) start(x))), seed.end = unlist(lapply(seed.match, function(x) end(x))))
	
	# Finally, write out the data frame
	out.name <- paste('Insert file path1', miRNA.list$names[miRNA], 'Insert file path2', sep = '')
	write.csv(full.match.df, out.name, row.names = FALSE, quote = FALSE)
}