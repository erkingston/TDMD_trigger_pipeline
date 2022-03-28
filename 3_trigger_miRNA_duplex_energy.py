import sys
import RNA
import pandas as pd

##############
# This python script takes as input the output file from 2_TDMD_trigger_narrower_looplength.R
# And it requires a path name for where the duplex energy data should be written out

infilename = sys.argv[1]
outfilename = sys.argv[2]

# Initialize all of the various parameters; this code pulls out the energy and pairing diagram for a given sequence
LastLine = False
sequences = open(infilename, 'rU')
energies = []
names = []
structures = []
seqs = []
header = sequences.readline()

# Manually specify the miRNA sequence for the miRNA of interest 
# Note 1: This should be written as DNA (so with Ts instead of Us)
# Note 2: Nucleotides 9-11 are represented by Ns, because we don't want to count pairing to those nucleotides towards the duplex energy of the site
DNAmiRNA = 'TATTGCACNNNTCCCGGCCTTT'
while LastLine == False:
	line = sequences.readline()
	if line == '':
		LastLine = True
	else:
		line = line.split(',')

		# Convert both the miRNA sequence (specified above) and the trigger sequence (from the file) to RNA
		miRNA = DNAmiRNA.replace('t','u')
		RNAline = line[2].replace('t', 'u')
	
		# Use duplexfold to calculate the energy of the trigger/miRNA duo, and also the pairing structure
		data = RNA.duplexfold(RNAline, miRNA)
		energy = data.energy
		energies.append(energy)
		names.append(line[0])
		structures.append(data.structure)
		seqs.append(line[2])

# Write all of the data out to the specified output file
mergeddata = list(zip(names, energies, structures, seqs))
df = pd.DataFrame(data = mergeddata, columns = ['names', 'energy', 'structure', 'sequence'])
df.to_csv(outfilename, index = False, header = False)
sequences.close()