'''
codonTable.py
'''

codonTable = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"X", "UAG":"X",
    "UGU":"C", "UGC":"C", "UGA":"X", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

newCodonTable = {}

for i in codonTable:
	newCodonTable[i] = codonTable[i]
	newCodonTable[i.replace('U','T')]=codonTable[i]

for i in newCodonTable.keys():
	newCodonTable[i[:2]+'N']='Z'
	newCodonTable[i[0]+'N'+i[2]]='Z'
	newCodonTable['N'+i[1:]]='Z'
	newCodonTable[i[0]+'NN']='Z'
	newCodonTable['N'+i[1]+'N']='Z'
	newCodonTable['NN'+i[2]]='Z'
	newCodonTable['NNN']='Z'

codonTable = newCodonTable
