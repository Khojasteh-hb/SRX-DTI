
import re
import readFasta
import numpy as np
import pandas as pd

def GTPC(fastas, **kw):
	group = {
		'alphaticr': 'GAVLMI',
		'aromatic': 'FYW',
		'postivecharger': 'KRH',
		'negativecharger': 'DE',
		'uncharger': 'STCPNQ'
	}

	groupKey = group.keys()
	baseNum = len(groupKey)
	triple = [g1+'.'+g2+'.'+g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]

	index = {}
	for key in groupKey:
		for aa in group[key]:
			index[aa] = key

	encodings = []
	header = ['#'] + triple
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])

		code = [name]
		myDict = {}
		for t in triple:
			myDict[t] = 0

		sum = 0
		for j in range(len(sequence) - 3 + 1):
			myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] = myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] + 1
			sum = sum +1

		if sum == 0:
			for t in triple:
				code.append(0)
		else:
			for t in triple:
				code.append(myDict[t]/sum)
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

result = GTPC(fastas1, **kw)

data1 = np.matrix(result[1:])[:, 1:]
data_gtpc = pd.DataFrame(data=data1)

data_gtpc.to_csv('EN_gtpc.csv')
#data_gtpc.to_csv('GPCR_gtpc.csv')
#data_gtpc.to_csv('IC_gtpc.csv')
#data_gtpc.to_csv('NR_gtpc.csv')