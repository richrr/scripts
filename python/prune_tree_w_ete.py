import os
import sys
from utils import *
from ete3 import Tree

# usage: cd /nfs3/PHARM/Morgun_Lab/richrr/ebiomedicine_review_paper/etetoolkit
# python ~/Morgun_Lab/richrr/scripts/python/prune_tree_w_ete.py Bacteria.nw Bacteria_pruned_no_species.nw
# python ~/Morgun_Lab/richrr/scripts/python/prune_tree_w_ete.py Bacteria.nw Bacteria_pruned_no_genus_species.nw


# http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees

# Load a tree structure from a newick file.
t = Tree(sys.argv[1], format=1)


P = list() # for pruning
#for level in ["phylum", "class", "order", "family", "genus"]: # "species"
for level in ["phylum", "class", "order", "family"]: # "genus", "species"
    print level
    tmp = t.search_nodes(rank=level)
    print len(tmp)
    tmp = [o.name for o in tmp]
    P.extend(tmp)

print "Nodes for pruning"
P = list(set(P))
print len(P)

#writeLIST_to_file(P, "nodes_Sent_for_pruning.txt")
writeLIST_to_file(P, "nodes_Sent_for_pruning_P_F.txt")  # including family. nothing from genus onwards

#phyl =  t.search_nodes(rank="phylum")
#t.prune([p.name for p in phyl])

'''
spec =  t.search_nodes(rank="species")
# use to remove_child()
'''

t.prune(P)
t.write(outfile="frmt0_"+sys.argv[2])
t.write(format=1, outfile="frmt1_"+sys.argv[2])


print "Done"

