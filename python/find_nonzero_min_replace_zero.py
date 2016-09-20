import os
import sys
from utils import *
import numpy as np

infile = sys.argv[1]
f = read_file(infile)

minn = 1
for i in f:
  for x in i.split('\t'): # each element in the line
    try:
        numb = float(x)   # try converting to float
        if numb < minn and numb != 0.0: # find min. non-zero float
            minn = numb
    except ValueError:
        continue

minn = minn/float(10)
print "Replacing 0.0 with %f\n" %minn

out_f = []
for i in f:
  temp_l = [str(minn) if x=='0.0' else x for x in i.strip('\n').split('\t')]
  out_f.append('\t'.join(temp_l))

outfile = infile.replace('.txt', '_non_zero.txt')
writeLIST_to_file(out_f, outfile)


