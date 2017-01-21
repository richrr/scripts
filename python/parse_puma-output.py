
'''
the taxa are compiled as per the occurances of |
since the data is relativized, we expect for each
taxa level, the sum of all taxa in each sample will be 100%
'''

"""
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
"""

import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re
import numpy as np

# Usage: python ~/scripts/python/collapse_same_taxons_into_single_row.py -i infile.txt -o summarized_taxon_out.txt -f 1 -t 9 -s 0


def main(args):

    parser = argparse.ArgumentParser(description='Summarize number of seqs in a sample having same taxonomy/OTU')
    parser.add_argument('-i', '--infile')
    parser.add_argument('-o', '--outfile', default="summarized_taxon_out.txt")  # output filename
    parser.add_argument('-f', '--fromcolumn', default=0, type=int) # 0-index based column from where the samples start
    parser.add_argument('-t', '--tillcolumn', default=99999, type=int) # 0-index based column where the samples end (inclusive)
    parser.add_argument('-s', '--summarizecolumn', nargs='*', type =int, default=[0])  # creates a list of items (0-index based column) used for summarizing
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    
    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile  = args.infile
    outfile = args.outfile
    fromcol = args.fromcolumn
    tillcol = args.tillcolumn
    summarcol = args.summarizecolumn
    delim = args.delimiter

    d_list=list()
    comment_list=list()

    for l in open(infile):
        if '#' not in l:
            d_list.append(l.strip('\n'))
        else:
            comment_list.append(l.strip())


    # since columns are skipped (#OTU ID, a string column) writeLIST_to_file(comment_list, outfile) will not work.
    cmt_temp_dict=dict()
    for l in comment_list:
        contents = l.strip('\n').split(delim)
        key = delim.join([contents[i] for i in summarcol])
        try:
            value= [contents[col] for col in np.arange(fromcol, tillcol+1)]
        except:
            print l
        if key not in cmt_temp_dict:
            cmt_temp_dict[key] = value
        else:
            print l
    out = open(outfile, 'w')
    for k,v in cmt_temp_dict.items():
        val = [str(j) for j in v]
        out.write ('%s%s%s\n' %(k, delim, delim.join(val)))
    out.close()


    temp_dict=dict()
    for l in d_list:
        if "root;" not in l:
            continue
        contents = l.strip('\n').split(delim)
        #print contents
        key = delim.join([contents[i] for i in summarcol])
        # replace any ':' in key with ';' to avoid where files came from edited biom files. it is optional but i recommend
        #key = key.replace(':', ';')
        key = key.replace('"' , '')
        #key = key.replace(' ', '') # remove spaces since the join cmd used downstream cannot insert tab and so I need another script to convert space to tab

        # use the second last element as "species"
        new_key = key.split(';')
        key = new_key[-2]
        
        try:
            value= [float(contents[col]) for col in np.arange(fromcol, tillcol+1)]
        except:
            print l
        if key not in temp_dict:
            temp_dict[key] = value
        else:
            temp_dict[key] = [x+y for x,y in zip(temp_dict[key], value)]
    # this might be a list of tuples [(k1,v1), (k2,v2)], where each v is a list [vi vii ...]
    sorted_temp_dict = sorted(temp_dict.iteritems(), key=operator.itemgetter(1), reverse=True) 
    out = open(outfile, 'a')
    for i in sorted_temp_dict:
        val = [str(j) for j in i[1]]
        out.write ('%s%s%s\n' %(i[0], delim, delim.join(val)))
    out.close()
    #print len(sorted_temp_dict)




if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

