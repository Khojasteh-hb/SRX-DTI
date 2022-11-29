
import re
import readFasta
import numpy as np
import pandas as pd

def DPC(fastas, **kw):
	AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWY'
	encodings = []
	diPeptides = [aa1 + aa2 for aa1 in AA for aa2 in AA]
	header = ['#'] + diPeptides
	encodings.append(header)

	AADict = {}
	for i in range(len(AA)):
		AADict[AA[i]] = i

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		tmpCode = [0] * 400
		for j in range(len(sequence) - 2 + 1):
			tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+1]]] = tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+1]]] +1
		if sum(tmpCode) != 0:
			tmpCode = [i/sum(tmpCode) for i in tmpCode]
		code = code + tmpCode
		encodings.append(code)
	return encodings

kw = {'path': r"Enzyme.txt", 'order': 'ACDEFGHIKLMNPQRSTVWY'}
#kw = {'path': r"GPCR.txt", 'order': 'ACDEFGHIKLMNPQRSTVWY'}
#kw = {'path': r"Ion channel.txt", 'order': 'ACDEFGHIKLMNPQRSTVWY'}
#kw = {'path': r"Nuclear receptor.txt", 'order': 'ACDEFGHIKLMNPQRSTVWY'}

fastas1 = readFasta.readFasta(r"Enzyme.txt")
#fastas1 = readFasta.readFasta(r"GPCR.txt")
#fastas1 = readFasta.readFasta(r"Ion channel.txt")
#fastas1 = readFasta.readFasta(r"Nuclear receptor.txt")

result = DPC(fastas1, **kw)

data1 = np.matrix(result[1:])[:, 1:]
data_dpc = pd.DataFrame(data=data1)

data_dpc.to_csv('EN_dpc.csv')
#data_dpc.to_csv('GPCR_dpc.csv')
#data_dpc.to_csv('IC_dpc.csv')
#data_dpc.to_csv('NR_dpc.csv')