import os
import sys
import string
import re
from utils import *

infile = sys.argv[1]
mapfile = sys.argv[2]

delim = '\t'
maplist = list_to_dict_all_vals(read_file(mapfile), delim, ';')

outputlist = list()
with open(infile) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if '#' in line:
            outputlist.append(line)
            continue
        cols = line.split(delim)
        new_line = line
        if cols[0] in maplist:
            new_line = "%s\t%s" %(line, maplist[cols[0]])
        outputlist.append(new_line)

outfile = "merged-" + infile.replace('../', '')
writeLIST_to_file(outputlist, outfile)

  

