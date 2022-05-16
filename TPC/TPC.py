
import re
import readFasta
import numpy as np
import pandas as pd

def TPC(fastas, **kw):
	AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWY'
	encodings = []
	triPeptides = [aa1 + aa2 + aa3 for aa1 in AA for aa2 in AA for aa3 in AA]
	header = ['#'] + triPeptides
	encodings.append(header)

	AADict = {}
	for i in range(len(AA)):
		AADict[AA[i]] = i

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		tmpCode = [0] * 8000
		for j in range(len(sequence) - 3 + 1):
			tmpCode[AADict[sequence[j]] * 400 + AADict[sequence[j+1]]*20 + AADict[sequence[j+2]]] = tmpCode[AADict[sequence[j]] * 400 + AADict[sequence[j+1]]*20 + AADict[sequence[j+2]]] +1
		if sum(tmpCode) != 0:
			tmpCode = [i/sum(tmpCode) for i in tmpCode]
		code = code + tmpCode
		encodings.append(code)
	return encodings

kw = {'path': r"Enzyme.txt",'order': 'ACDEFGHIKLMNPQRSTVWY'}
fastas1 = readFasta.readFasta(r"Enzyme.txt")

result = TPC(fastas1, **kw)

data1 = np.matrix(result[1:])[:,1:]
data_ = pd.DataFrame(data=data1)
data_.to_csv('TPC_Enzyme.csv')