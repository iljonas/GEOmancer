import csv, os, gzip
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
masterPath = "C:\\Users\\" + os.getlogin() + "\\Documents\\Capstone Files\\GEO-Antimicrobial-Adjunct-Project\\Output_Files\\Exports\\output.tsv.gz"

with gzip.open(masterPath, 'rt') as mFile:
	mReader = mFile.readlines()
	masterData = [[col for col in row.rstrip('\n').split('\t')] for row in mReader]
mFile.close()

for mRow in masterData:
	if mRow[1] == '"Gene"':
		geneIDs.append(mRow[0].replace('"', ''))
	else:
		proteinIDs.append(mRow[0].replace('"', ''))

cols = 'id,reviewed,protein%20names,families,organism,organism-id,sequence,genes'

cLen = cols.count(',') + 1

if(len(geneIDs) > 0):
	geneIDs = ' '.join(geneIDs)
	geneData = list(csv.reader(getUniProtData(geneIDs, cols, 'GENENAME'), delimiter = '\t'))
	newGene = ['\t'.join(geneRow[-2] + geneRow[0:-3]) for geneRow in geneData if len(geneRow) > 2]
	print(newGene)

# if(len(proteinIDs) > 0):
	# proteinIDs = ' '.join(proteinIDs)
	# protData = list(csv.reader(getUniProtData(proteinIDs, cols, 'ACC+ID'), delimiter = '\t'))
	# newProt = [protRow[-2] + protRow[0:cLen] for protRow in protData]
	# obsIDs = []

	# for i in range(len(protData) - 1):
		# if not protData[i][1]:
			# print(protData[i][2])
			# protData[i][2] = protData[i][2][len("Merged into ")  : -1]
			# obsIDs.append(protData[i][0,2])
			# del protData[i]

	# updatedData = list(csv.reader(getUniProtData(obsIDs[1], cols, 'ACC+ID'), delimiter = '\t'))
	# newUpdate = []
	# for obsRow in obsIDs:
		# for updatedRow in updatedData:
			# if obsIDs[1] == updatedRow[0]:
				# newUpdate.append(obsIDs[0] + updatedRow)
				# break

newGene[0][0] = "Original Entry"
with open("C:\\Users\\iljonas\\Desktop\\testFile.txt", 'w') as outFile:
	outFile.write('\n'.join(['\t'.join([str(element) for element in row]) for row in geneData]))
#	outFile.write(newProt[1:])
#	outFile.write(newUpdate[1:])
outFile.close()
print(page)


#So um, query for genes to ACC and ACC+ID to ACC
#From that second list, look for any status entries that are blank, as this will indicate if they're obsolete
#	If the above isn't true, look for any protein names that start with "Merged into"
#Take the list ID from "yourlist:" and put it into the uniprot URL:
#	http://www.uniprot.org/uniprot/?query=yourlist:M20180621A7434721E10EE6586998A056CCD0537E7F9966Y&format=tab
#Rerun ACC+ID query for the updated proteins under protein names
#To get all proteins:
#https://www.uniprot.org/uniprot/?query=yourlist:M20180621A7434721E10EE6586998A056CCD0537E7F9966Y&columns=id,reviewed,protein%20names,families,organism,organism-id,sequence,genes&format=tab