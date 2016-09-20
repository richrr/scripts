import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re
import numpy as np

# Usage: python ~/scripts/python/summarize_topk_taxon.py -i top-10k-otus-taxons-rosana-16srrna-row_sum_sanitized_summary_otu_table_mc2_w_tax_no_pynast_failures.biom.txt-rep_set_tax_assign.txt -o summarize_topk_taxon_out.txt -f 5 -t 7 -s 0 1

########## add code to print header row in the output file

def main(args):

    parser = argparse.ArgumentParser(description='Summarize number of seqs in a sample having same taxonomy/OTU')
    parser.add_argument('-i', '--infile')
    parser.add_argument('-o', '--outfile', default="summarize_topk_taxon_out.txt")  # output filename
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
    temp_dict=dict()

 
    for l in open(infile):
        if '#' not in l:
            d_list.append(l.strip('\n'))
        else:
            comment_list.append(l.strip())

    '''
    for col in np.arange(fromcol, tillcol+1):
        temp_dict=dict()
        for l in d_list:
            contents = l.strip().split(delim)
            key = delim.join([contents[i] for i in summarcol])
            #key = key.replace('[' , '_').replace(']', '_')
            if key not in temp_dict:
                temp_dict[key] = contents[col]
            else:
                temp_dict[key] = temp_dict[key] + contents[col]
            # list of tuples, ordered as per the index/value of the dictionary
        sorted_temp_dict = sorted(temp_dict.iteritems(), key=operator.itemgetter(1))
        print [i for i in sorted_temp_dict[:5]]
    '''

    for l in d_list:
        contents = l.strip('\n').split(delim)
        #print contents
        key = delim.join([contents[i] for i in summarcol])
        try:
            value= [float(contents[col]) for col in np.arange(fromcol, tillcol+1)]
        except:
            print l
        if key not in temp_dict:
            temp_dict[key] = value
        else:
            temp_dict[key] = [x+y for x,y in zip(temp_dict[key], value)]
    sorted_temp_dict = sorted(temp_dict.iteritems(), key=operator.itemgetter(1), reverse=True)
    out = open(outfile, 'w')
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

