# adapted from uparse's uc2otutab.py
import os
import sys

# Create OTU table: 
#combine the otu-mappings from the previous step into an otu-by-sample table.
#This relies on the sequence ID incorporating the sample ID, as in 2_3456q being sample 2.
#This create a table with sample IDs for column headers, OTU IDS for row names, and counts in each cell.

# usage: python ~/scripts/python/uc_to_otu_table.py map.uc > otu_table.txt

FileName = sys.argv[1]

OTUIds = []
SampleIds = []
OTUTable = {}


#http://drive5.com/usearch/manual/ucout.html
with open(FileName) as f:
    for l in f:
        l = l.strip()
        if l[0] != 'H':
            continue
        recs = l.split('\t')
        querylabel = recs[8]
        OTUId = recs[9]
        if OTUId not in OTUIds:
	    OTUIds.append(OTUId)
	    OTUTable[OTUId] = {}
        SampleId = querylabel[:querylabel.find('_')]
        if SampleId not in SampleIds:
	    SampleIds.append(SampleId)
	try:
	    OTUTable[OTUId][SampleId] += 1
	except:
	    OTUTable[OTUId][SampleId] = 1
f.close()

s = "OTUId"
for SampleId in SampleIds:
	s += "\t" + SampleId
print s

for OTUId in OTUIds:
	s = OTUId
	for SampleId in SampleIds:
		try:
			n = OTUTable[OTUId][SampleId]
		except:
			n = 0
		s += "\t" + str(n)
	print s
