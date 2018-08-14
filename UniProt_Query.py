import csv, os, gzip
from sys import argv
from urllib.parse import urlencode
from urllib.request import Request, urlopen, getproxies

def getUniProtData(query, cols, fType, tType = 'ACC'):
	url = 'https://www.uniprot.org/uploadlists/'
	params = {
		'from':fType,
		'to':tType,
		'active':'no',
		'format':'tab',
		'query':query,
		'columns':cols
	}
	data = urlencode(params).encode("utf-8")
	req = Request(url, data)
	contact = "iljonas6@gmail.com" # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
	req.add_header('User-Agent', 'Python %s' % contact)
	response = urlopen(req)
	return(response.read().decode('utf-8').split('\n'))

geneIDs = []
proteinIDs = []
masterPath = ['C:', 'Users', os.getlogin(), 'Documents', 'Capstone Files', 'GEO-Antimicrobial-Adjunct-Project', 'Output_Files', 'Exports']
masterFile = os.sep.join(masterPath) + os.sep + argv[-1] + '.tsv.gz'

with gzip.open(masterFile, 'rt') as mFile:
	mReader = mFile.readlines()
	masterData = [[col for col in row.rstrip('\n').split('\t')] for row in mReader]
mFile.close()

for mRow in masterData[1:]:
	if mRow[1] == '"Gene"':
		geneIDs.append(mRow[0].replace('"', ''))
	elif mRow[1] == '"Protein"':
		proteinIDs.append(mRow[0].replace('"', ''))

cols = 'id,reviewed,protein%20names,families,organism,organism-id,sequence,genes'
cLen = cols.count(',') + 1

if(len(geneIDs) > 0):
	geneIDs = ' '.join(geneIDs)
	print('Querying UniProt for project genes...')
	geneData = list(csv.reader(getUniProtData(geneIDs, cols, 'GENENAME'), delimiter = '\t'))
	uniData = [(geneRow[-2] + '\t' + '\t'.join(geneRow[0:-2])) for geneRow in geneData if len(geneRow) > 2]

if(len(proteinIDs) > 0):
	proteinIDs = ' '.join(proteinIDs)
	print('Querying UniProt for project proteins...')
	protData = list(csv.reader(getUniProtData(proteinIDs, cols, 'ACC+ID'), delimiter = '\t'))
	obsIDs = []

	for i in range(len(protData) - 1):
		if not protData[i][1]:
			print(protData[i][2])
			protData[i][2] = protData[i][2][len("Merged into ")  : -1]
			obsIDs.append(protData[i][0,2])
			del protData[i]
	
	obsIDs_string = ' '.join([obsRow[1] for obsRow in obsIDs])
	print('Resubmitting outdated protein IDs to UniProt with their current ID...')
	updatedData = list(csv.reader(getUniProtData(obsIDs_string, cols, 'ACC+ID'), delimiter = '\t'))
	#newUpdate = []
	for obsRow in obsIDs:
		for updatedRow in updatedData:
			if obsIDs[1] == updatedRow[0]:
				protData.extend(obsIDs[0] + updatedRow[1:])
				#newUpdate.append(obsIDs[0] + updatedRow)
				break

	#newProt = [(protRow[-2] + '\t' + '\t'.join(protRow[:-2])) for protRow in protData if len(protRow) > 2]
	uniData.extend([(protRow[-2] + '\t' + '\t'.join(protRow[:-2])) for protRow in protData[1:] if len(protRow) > 2])
	
print('UniProt queries completed. Saving results')
uniFile = os.sep.join(masterPath[:-1]) + os.sep + 'UniProt_Downloads' + os.sep + 'uniProt-' + argv[-1] + '.tsv.gz'
with gzip.open(uniFile, 'wt') as outFile:
	outFile.write(uniData[0].replace(uniData[0][0:uniData[0].find('\t')], 'Expression ID') + '\n')
	outFile.write('\n'.join(uniData[1:]))
outFile.close()