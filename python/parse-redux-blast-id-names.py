import os
import sys
import string
import re
from utils import *

infile = sys.argv[1]

query = otu_id = read_info = ''
counter = 0
taxa_list = list()
otu_id_taxa_dict = dict()
#end_flag = False
with open(infile) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if "Query=" in line:
            query, otu_id, read_info = line.split(' ')
            counter = 0
            continue
        if '--' in line:
            if otu_id in otu_id_taxa_dict:
                sys.exit("this otu id already exists!")
            otu_id_taxa_dict[otu_id] = '\t'.join(taxa_list)
            query = otu_id = read_info = ''
            counter = 0
            taxa_list = list()
            continue
        if "ref|" in line and counter < 3:
            #ref_name_list = line.split('  ')
            line = line.replace('  ', ' ')
            taxa_list.append(line)
            counter += 1
            continue


if otu_id in otu_id_taxa_dict:
    sys.exit("this otu id already exists!")

otu_id_taxa_dict[otu_id] = '\t'.join(taxa_list)

#print otu_id_taxa_dict
writeDICT_to_file(otu_id_taxa_dict, "output-"+infile)


'''
Query= 104129 IR10_S_R12__68184

Length=226


     
Sequences producing significant alignments:  

ref|NZ_AXAI01000013.1|  Bradyrhizobium sp. Tv2a-2 A3AIDRAFT_sc...
ref|NZ_AWZU01000001.1|  Bradyrhizobium sp. ARR65 BraARR65DRAFT...
ref|NZ_AXAU01000005.1|  Bradyrhizobium elkanii WSM1741 YUODRAF...
ref|NZ_AXAP01000139.1|  Bradyrhizobium elkanii WSM2783 YY7DRAF...
ref|NZ_AXAP01000127.1|  Bradyrhizobium elkanii WSM2783 YY7DRAF...
ref|NZ_AXAP01000126.1|  Bradyrhizobium elkanii WSM2783 YY7DRAF...
ref|NZ_AXAP01000020.1|  Bradyrhizobium elkanii WSM2783 YY7DRAF...
ref|NZ_JNIJ01000042.1|  Bradyrhizobium sp. URHD0069 N554DRAFT_...
--
Query= 107036 IR10_S_R12__17110


'''
