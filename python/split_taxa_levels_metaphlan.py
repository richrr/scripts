import os
import sys
from utils import *

'''
the taxa are compiled as per the occurances of |
since the data is relativized, we expect for each
taxa level, the sum of all taxa in each sample will be 100%
'''
def check_taxa_level(l):
    flag = True
    #l_first = 
    numb_samples = len(l[0].split('\t'))
    if numb_samples == 0:
        sys.exit("Error in the list being passed")

    sum_list = [0 for i in range(numb_samples)]

    for row in l:
        cont = row.split('\t')
        if len(cont) != numb_samples:
            continue
        #print cont
        l_row = [float(i) for i in cont[1:]]
        sum_list = map(sum,zip(sum_list, l_row))
    
    #print sum_list
    
    for val in sum_list:
        if int(round(val)) != 100:
            flag = False
    return flag

infile = sys.argv[1]

headers = list()

taxa_levels = [["#k"], ["#p"],["#c"],["#o"],["#f"],["#g"],["#s"], ["#t"]]

with open(infile) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if '#' in line:
            headers.append(line)
            continue
        contents = line.split('\t')
        indx = len(contents[0].split('|')) - 1
        taxa_levels[indx].append(line.replace('|' , ';'))
        #sys.exit(0)


for l in taxa_levels:
    #print l
    if len(l) == 1: # there is only 1 element (e.g. "#t") in this list
        continue
    flag = check_taxa_level(l[1:]) # skip the taxa level indicator
    if (flag):
        new_l = [l[0]] + headers + l[1:]
        writeLIST_to_file(new_l, infile+"__"+l[0].replace('#','')+".tsv")
    else:
        new_l = [l[0]] + headers + l[1:]
        writeLIST_to_file(new_l, infile+"__"+l[0].replace('#','')+"-doesnt-sum-to-100.tsv")        
        sys.exit("Error: the values for %s per sample do not sum to 100" %l[0])

