import os
import sys
from utils import *

# usage: 
# cd /nfs3/PHARM/Morgun_Lab/richrr/ebiomedicine_review_paper/ncbi_taxonomy_new_dump_Mar_6_2019
# python ~/Morgun_Lab/richrr/scripts/python/ncbi_nodes_to_names.py levels_to_keep.txt scientific_names.dmp keep_names_for_pruning.txt

cont = read_file(sys.argv[1])
vals = list()
for l in cont:
    l = l.strip('\n')
    if l: # accounts for empty line
        #print l
        c = l.split('|')
        vals.append(c[0].strip())

print len(vals)

names = read_file(sys.argv[2])
id_name_dict = dict()
for n in names:
    n = n.strip('\n')
    if 'scientific name' in n:
        c = n.split('\t|\t')
        id_name_dict[c[0]] = c[1]

print len(id_name_dict.keys())

outids = list()
for v in vals:
    outids.append(id_name_dict[v])

writeLIST_to_file(outids, sys.argv[3])

print "Done"