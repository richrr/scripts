import os
import sys
from utils import *

# this code takes a node file and makes a new node for each type mentioned and makes
# the remaining nodes in that type to ignore. mostly useful to make per pheno node files

### see if/how this code is different than the make-bibc-cmds-per-node.R ###

node_file = sys.argv[1]
out_node_file = node_file.split('/')[-1]
print node_file

type = sys.argv[2]
print type

edge_file = sys.argv[3]
print edge_file

lines = read_file(node_file)

as_is = list() # these lines will not be changed in the output files
change_nodes = list() # each of these nodes will have a file in which they have the type and ignore in other files

for l in lines:
    l = l.strip()
    if not l:
        continue
    conts = l.split(',')
    node = conts[0]
    categ = conts[1]
    #print node
    #print categ
    if categ == type:
        change_nodes.append(node)
    else:
        as_is.append(l)

cmds_list = list()

for n in change_nodes:
    print n
    tmp_list=list()
    tmp_list.extend(as_is)
    tmp_list.append(n + ',' + type)
    rem_change_nodes = [ c + ',ignore' for c in change_nodes if c != n ]
    #print rem_change_nodes
    tmp_list.extend(rem_change_nodes)
    n_file_name = ''.join(e for e in n if e.isalnum())
    new_file_name = out_node_file + '.' + n_file_name
    
    writeLIST_to_file(tmp_list, new_file_name)
    
    cmds1 = ' '.join(['Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/prep-for-betw-central.R', edge_file, new_file_name, new_file_name])
    cmds2 = 'SGE_Batch -c "/nfs3/PHARM/Morgun_Lab/richrr/scripts/R/betweennessadapter node_type microbe degs ./' + new_file_name + '.graphml > bibc-' + new_file_name + '-results.txt" -m 50G -F 50G -q samwise -M rodrrich@oregonstate.edu -r log_' + new_file_name
    cmds_list.append(cmds1)
    cmds_list.append(cmds2)
    

writeLIST_to_file(cmds_list, out_node_file + '-cmds.txt')


