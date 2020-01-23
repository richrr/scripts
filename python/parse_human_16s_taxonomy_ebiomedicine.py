import os
import sys
from utils import *

# cd /nfs3/PHARM/Morgun_Lab/richrr/ebiomedicine_review_paper/human_closed_ref_t2d_16s_data
# python ~/Morgun_Lab/richrr/scripts/python/parse_human_16s_taxonomy_ebiomedicine.py uniq_otus_detected_in_humans.txt /nfs3/PHARM/Morgun_Lab/richrr/db/amiran_mothur_format/qiime_otu_top_blast.csv

cont = read_file(sys.argv[1])
keep = [c.strip('\n') for c in cont]
print len(keep)


names = read_file(sys.argv[2])
id_name_dict = list_to_dict(names, delim=',', joiner=',', string="current", key_column=0, val_column=1)
print len(id_name_dict.keys())

dictkeys = id_name_dict.keys()
# remove otus not in human data
for key in dictkeys:
    if key not in keep:
        del id_name_dict[key]

print len(id_name_dict.keys())
taxa = set(id_name_dict.values())
print len(taxa)

writeDICT_to_file(id_name_dict, sys.argv[1]+".w.taxa.tsv", delim='\t')


print "Done"
