import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re
import numpy as np
import matplotlib
matplotlib.use('Agg') # to avoid server trying to display the figures
import matplotlib.pyplot as plt

# Usage: python ~/scripts/python/summarize_topk_taxon.py -i top-10k-otus-taxons-rosana-16srrna-row_sum_sanitized_summary_otu_table_mc2_w_tax_no_pynast_failures.biom.txt-rep_set_tax_assign.txt -o summarize_topk_taxon_out.txt -f 5 -t 7 -s 0 1

########## add code to print header row in the output file

def main(args):

    parser = argparse.ArgumentParser(description='Walk through the sorted otu lists and plot number of overlapping same taxonomy/OTU')
    parser.add_argument('-i', '--infile1')
    parser.add_argument('-j', '--infile2')
    parser.add_argument('-o', '--outfile', default="out.jpg")  # output filename
    #parser.add_argument('-f', '--fromcolumn', default=0, type=int) # 0-index based column from where the samples start
    #parser.add_argument('-t', '--tillcolumn', default=99999, type=int) # 0-index based column where the samples end (inclusive)
    #parser.add_argument('-s', '--summarizecolumn', nargs='*', type =int, default=[0])  # creates a list of items (0-index based column) used for summarizing
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    
    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile1  = read_file(args.infile1)
    infile2  = read_file(args.infile2)
    outfile = args.outfile
    delim = args.delimiter



    otu_id_col_file1 = taxonomy_col_file1 = otu_id_col_file2 = taxonomy_col_file2 = '' # ids for the columns

    file1_headers_list = list()
    file2_headers_list = list()

    file1_dict = dict()
    file2_dict = dict()

    file1_otusum_dict = dict()
    file2_otusum_dict = dict()

    for l in infile1:
        l = l.strip('\n')
        if not l:
            continue
        contents = l.split(delim)
        if '#' in l:
            if '#OTU ID' in contents:
                otu_id_col_file1 = contents.index('#OTU ID')
            if 'taxonomy' in contents:
                taxonomy_col_file1 = contents.index('taxonomy')
            file1_headers_list.append(l)
            continue
        
        key = contents[otu_id_col_file1]
        value = ''
        try:
            value= [float(contents[col]) for col in np.arange(otu_id_col_file1+1, taxonomy_col_file1)] # start after otu id and end before taxonomy
        except:
            print l
            sys.exit("A")
        if key not in file1_dict:
            file1_otusum_dict[key] = sum(value)
            value.append(contents[taxonomy_col_file1]) # add taxonomy to the list
            file1_dict[key] = value
        else:
            print l
            sys.exit("B")
        

    for l in infile2:
        l = l.strip('\n')
        if not l:
            continue
        contents = l.split(delim)
        if '#' in l:
            if '#OTU ID' in contents:
                otu_id_col_file2 = contents.index('#OTU ID')
            if 'taxonomy' in contents:
                taxonomy_col_file2 = contents.index('taxonomy')
            file2_headers_list.append(l)
            continue
        
        key = contents[otu_id_col_file2]
        value = ''
        try:
            value= [float(contents[col]) for col in np.arange(otu_id_col_file2+1, taxonomy_col_file2)] # start after otu id and end before taxonomy
        except:
            print l
            sys.exit("C")
        if key not in file2_dict:
            file2_otusum_dict[key] = sum(value)
            value.append(contents[taxonomy_col_file2]) # add taxonomy to the list
            file2_dict[key] = value
        else:
            print l
            sys.exit("D")       

    # sorts dict as per value, returns a list of tuples (key, value)
    # for cases where values are same for multiple keys, the keys stay the same as in their original order (not the best)
    #sorted_temp_dict1 = sorted(file1_otusum_dict.iteritems(), key=operator.itemgetter(1), reverse=True)
    #sorted_temp_dict2 = sorted(file2_otusum_dict.iteritems(), key=operator.itemgetter(1), reverse=True)

    #http://stackoverflow.com/questions/33195384/python-sorting-a-dictionary-based-on-values-and-order-based-on-key-if-values-ar

    ListA = list()
    ListB = list()

    out = open(args.infile1+".sorted.txt", 'w')
    for l in file1_headers_list:
        out.write ('%s\n' %(l))

    #for i in sorted_temp_dict1:
    #    k = i[0]
    # for cases where values are same for multiple keys, the keys are asc sorted in lexical order
    for k in sorted(file1_otusum_dict, key=lambda x: (-file1_otusum_dict.get(x), x)):
        ListA.append(k)
        val = [str(j) for j in file1_dict[k]]
        out.write ('%s%s%s\n' %(k, delim, delim.join(val)))
    out.close()

    out = open(args.infile2+".sorted.txt", 'w')
    for l in file2_headers_list:
        out.write ('%s\n' %(l))

    #for i in sorted_temp_dict2:
    #    k = i[0]
    # for cases where values are same for multiple keys, the keys are asc sorted in lexical order
    for k in sorted(file2_otusum_dict, key=lambda x: (-file2_otusum_dict.get(x), x)):
        ListB.append(k)
        val = [str(j) for j in file2_dict[k]]
        out.write ('%s%s%s\n' %(k, delim, delim.join(val)))
    out.close()

    #JIvalues is a list of pairs [x,y], where 'y' is the JI at cutoff 'x'
    JIvalues = jaccardIndex(ListA, ListB)
    plt.ioff()
    '''
    plt.scatter(*zip(*JIvalues))
    plt.savefig('miseq-JI.png')
    plt.close()
    '''
    plt.figure(1)
    plt.subplot(211)
    plt.scatter(*zip(*JIvalues))
    plt.ylabel('Jaccard Index')

    JIvalues = jaccardIndex_counts(ListA, ListB)
    '''
    plt.ioff()
    plt.scatter(*zip(*JIvalues))
    plt.savefig('miseq-count.png')
    plt.close()
    '''
    plt.subplot(212)
    plt.scatter(*zip(*JIvalues))
    plt.xlabel('Steps')
    plt.ylabel('Overlap')
    plt.savefig('miseq-JI_count_combo.png')
    plt.close()




if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

