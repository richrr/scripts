import sys
import os
from utils import *
from collections import OrderedDict

infile = sys.argv[1]
node_module_dict1 = readDict(infile, 1, 7, '\t')

res = OrderedDict()
for v, k in node_module_dict1.items(): # k is module, v is OTU
    if k in res:
        res[k].append(v)
    else:
        res[k] = [v]
writeDICT_to_file(res, sys.argv[2])


#http://stackoverflow.com/questions/6602881/text-file-to-list-in-r
