
import re
import readFasta
from collections import Counter
import numpy as np
import pandas as pd

def AAC(fastas, **kw):
	AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWY'
	#AA = 'ARNDCQEGHILKMFPSTWYV'
	encodings = []
	header = ['#']
	for i in AA:
		header.append(i)
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		count = Counter(sequence)
		for key in count:
			count[key] = count[key]/len(sequence)
		code = [name]
		for aa in AA:
			code.append(count[aa])
		encodings.append(code)
	return encodings

kw=  {'path': r"Enzyme.txt",'order': 'ACDEFGHIKLMNPQRSTVWY'}
fastas1 = readFasta.readFasta(r"Enzyme.txt")

result = AAC(fastas1, **kw)

data1 = np.matrix(result[1:])[:,1:]
data_ = pd.DataFrame(data=data1)
data_.to_csv('AAC_Enzyme.csv')