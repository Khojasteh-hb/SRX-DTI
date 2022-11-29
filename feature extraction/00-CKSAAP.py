
import sys, os
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import readFasta
import checkFasta
import numpy as np
import pandas as pd


USAGE = """
USAGE:
	python CKSAAP.py input.fasta <k_space> <output>

	input.fasta:      the input protein sequence file in fasta format.
	k_space:          the gap of two amino acids, integer, defaule: 5
	output:           the encoding file, default: 'encodings.tsv'
"""

def CKSAAP(fastas, gap=5, **kw):
	if gap < 0:
		print('Error: the gap should be equal or greater than zero' + '\n\n')
		return 0

	if checkFasta.minSequenceLength(fastas) < gap+2:
		print('Error: all the sequence length should be larger than the (gap value) + 2 = ' + str(gap+2) + '\n\n')
		return 0

	AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWY'
	encodings = []
	aaPairs = []
	for aa1 in AA:
		for aa2 in AA:
			aaPairs.append(aa1 + aa2)
	header = ['#']
	for g in range(gap+1):
		for aa in aaPairs:
			header.append(aa + '.gap' + str(g))
	encodings.append(header)
	for i in fastas:
		name, sequence = i[0], i[1]
		code = [name]
		for g in range(gap+1):
			myDict = {}
			for pair in aaPairs:
				myDict[pair] = 0
			sum = 0
			for index1 in range(len(sequence)):
				index2 = index1 + g + 1
				if index1 < len(sequence) and index2 < len(sequence) and sequence[index1] in AA and sequence[index2] in AA:
					myDict[sequence[index1] + sequence[index2]] = myDict[sequence[index1] + sequence[index2]] + 1
					sum = sum + 1
			for pair in aaPairs:
				code.append(myDict[pair] / sum)
		encodings.append(code)
	return encodings


#kw = {'path': r"Enzyme.txt", 'order': 'ACDEFGHIKLMNPQRSTVWY'}
#kw = {'path': r"GPCR.txt", 'order': 'ACDEFGHIKLMNPQRSTVWY'}
#kw = {'path': r"Ion channel.txt", 'order': 'ACDEFGHIKLMNPQRSTVWY'}
kw = {'path': r"Nuclear receptor.txt", 'order': 'ACDEFGHIKLMNPQRSTVWY'}

#fastas1 = readFasta.readFasta(r"Enzyme.txt")
#fastas1 = readFasta.readFasta(r"GPCR.txt")
#fastas1 = readFasta.readFasta(r"Ion channel.txt")
fastas1 = readFasta.readFasta(r"Nuclear receptor.txt")

result = CKSAAP(fastas1, **kw)

data1 = np.matrix(result[1:])[:,1:]
data_cksaap = pd.DataFrame(data=data1)

#data_cksaap.to_csv('EN_cksaap.csv')
#data_cksaap.to_csv('GPCR_cksaap.csv')
#data_cksaap.to_csv('IC_cksaap.csv')
data_cksaap.to_csv('NR_cksaap.csv')