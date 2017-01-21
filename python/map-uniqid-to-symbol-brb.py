import sys
import os
from utils import *


file1 = sys.argv[1]
file2 = sys.argv[2]
colf1 = int(sys.argv[3])
colf2 = int(sys.argv[4])
ofile = sys.argv[5]

ret_list = map_identifiers(file1, file2, "," , ",", colf1, colf2) 
writeLIST_to_file(ret_list, ofile)

