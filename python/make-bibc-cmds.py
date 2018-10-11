import os
import sys
from utils import *

# this code takes an edge file and a 2 column tab-delimited file consisting of "name" and "node file" per line
# and creates the bibc commands for them.

edge_file = sys.argv[1]
print edge_file


node_file = sys.argv[2]
node_dict = list_to_dict(read_file(node_file))
print node_dict

# to do from here #

cmds_list = list()

for k,v in node_dict.items():
    print k,v
    v = v.strip()

    
    cmds1 = ' '.join(['Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/prep-for-betw-central.R', edge_file, v, k])
    cmds2 = 'SGE_Batch -c "/nfs3/PHARM/Morgun_Lab/richrr/scripts/R/betweennessadapter node_type microbe degs ./' + k + '.graphml > bibc-' + k + '-results.txt" -m 50G -F 50G -q transkingdom -M rodrrich@oregonstate.edu -r log_' + k
    cmds_list.append(cmds1)
    cmds_list.append(cmds2)
    

writeLIST_to_file(cmds_list, 'bibc-cmds.txt')





