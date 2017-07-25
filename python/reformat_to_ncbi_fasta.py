import os
import sys
from utils import *

infile="/nfs3/PHARM/Morgun_Lab/richrr/db/amiran_mothur_format/ncbi20.fasta"
taxafile="/nfs3/PHARM/Morgun_Lab/richrr/db/amiran_mothur_format/ncbi20.tax"
n = 80

# create id to taxa mapping
taxainfo = read_file(taxafile)
id_taxa_dict = list_to_dict(taxainfo)


outfile = list()
seq=''
with open(infile) as f:
    for line in f:
        line=line.strip()
        #print line
        if line.startswith(">"):
            #print line
            # write the old sequence and resets the sequence
            strs = [seq[i:i+n] for i in range(0, len(seq), n)]
            # do not add empty lines
            if len(strs) > 0:
                outfile.extend(strs)
            newid = line.replace(">", "")
            # add taxa info
            if newid in id_taxa_dict:
                line = line + ' ' + id_taxa_dict[newid]
            outfile.append(line)
            print(line)
            seq=''
        else:
            seq = seq+line


writeLIST_to_file(outfile, infile+".ncbi.format.fa")
