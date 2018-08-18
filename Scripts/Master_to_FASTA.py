import os, gzip
from sys import argv
#from Bio.Blast import NCBIWWW

masterName = argv[-1]
if not masterName.lower().startswith('master-'):
	masterName = 'master-' + master

masterFolder = ['C:', 'Users', os.getlogin(), 'Documents', 'GEOmancer', 'Output_Files', 'Exports']
masterPath = os.sep.join(masterFolder) + os.sep + masterName + '.tsv.gz'

with gzip.open(masterPath, 'rt') as mFile:
	mReader = mFile.readlines()
	masterData = [[col for col in row.rstrip('\n').split('\t')] for row in mReader]
mFile.close()

geneIndex = masterData[0].index('"Expression.ID"')
seqIndex = masterData[0].index('"Sequence"')

fastaFile = ''
for row in masterData[1:]:
	fastaFile += '> ' + row[geneIndex].replace('"', '') + '\n' + row[seqIndex].replace('"', '') + '\n'
	
print(fastaFile)

fFile = os.sep.join(masterFolder[:-1]) + os.sep + 'UniProt_Downloads' + os.sep + argv[-1] + '.fasta'
with open(fFile, 'wt') as outFile:
	outFile.write(fastaFile)
outFile.close()

#resultHandle = NCBIWWW.qblast('blastp', 'nr', proteinList)
