import os, gzip, requests
from sys import argv
from urllib.parse import urlencode
from urllib.request import urlopen

def masterToFASTA(data, id)
	idIndex = data[0].index(id)
	seqIndex = data[0].index('"Sequence"')

	fasta = ''
	for row in data[1:]:
		if row[idIndex] != 'NA':
			fasta += '> ' + row[idIndex].replace('"', '') + '\n' + row[seqIndex].replace('"', '') + '\n'
	return(fasta)

def getCDHIT_Data(query, annotate = 0, cutOff = 0.9, globalSeq = 0, bestClust = 0, bandAlign = 20, \ 
	minAlign_long = 0.0, maxUnalign_long = 'unlimited', minAlign_short = 0.0, maxUnalign_short = 'unlimited', minSimilar = 0.0, maxDiff = 'unlimited'):

	url = 'http://weizhong-cluster.ucsd.edu/cdhit_suite/cgi-bin/index.cgi?cmd=cd-hit'
	files = {'SeqF':query}
	params = {
		#Incorporate annotation info at header line
		'anno':annotate
		#Sequence identity cut-off
		'lc1':cutOff
		#-G: use global sequence identity
		'uG':globalSeq
		#-g: sequence is clustered to the best cluster that meet the threshold
		'lg':bestClust
		#-b: bandwidth of alignment
		'lb':bandAlign
		#-aL: minimal alignment coverage (fraction) for the longer sequence
		'lauL':minAlign_long,
		#-AL: maximum  unaligned part (amino acids/bases) for the longer sequence
		'uAuL':maxUnalign_long,
		#-aS: minimal alignment coverage (fraction) for the shorter sequence
		'lauS':minAlign_short,
		#-AS: maximum  unaligned part (amino acids/bases) for the shorter sequence
		'uAuS':maxUnalign_short,
		#-s:  minimal length similarity (fraction)
		'ls':minSimilar,
		#-S:  maximum  length difference in amino acids/bases(-S)
		'uS':maxDiff
	}
	data = urlencode(params).encode("utf-8")
	req = requests.post(url, files=files, data=data)
	contact = "iljonas6@gmail.com" # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
	req.add_header('User-Agent', 'Python %s' % contact)
	response = urlopen(req)
	return(response.read().decode('utf-8').split('\n'))

masterName = argv[-1]
if not masterName.lower().startswith('master-'):
	masterName = 'master-' + master

masterFolder = ['C:', 'Users', os.getlogin(), 'Documents', 'GEOmancer', 'Output_Files', 'Exports']
masterPath = os.sep.join(masterFolder) + os.sep + masterName + '.tsv.gz'

with gzip.open(masterPath, 'rt') as mFile:
	mReader = mFile.readlines()
	masterData = [[col for col in row.rstrip('\n').split('\t')] for row in mReader]
mFile.close()

fastaFile = masterToFASTA(masterData, '"Expression.ID"')

fFile = os.sep.join(masterFolder[:-1]) + os.sep + 'UniProt_Downloads' + os.sep + argv[-1] + '.fasta'
with open(fFile, 'wt') as outFile:
	outFile.write(fastaFile)
outFile.close()