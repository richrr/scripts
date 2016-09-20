import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re

# Usage: python ~/scripts/python/topk_taxon.py -i row_sum_sanitized_summary_otu_table_mc2_w_tax_no_pynast_failures.biom.txt-rep_set_tax_assign.txt -o top-10k-otus-taxons-row_sum_sanitized_summary_otu_table_mc2_w_tax_no_pynast_failures.biom.txt-rep_set_tax_assign.txt -k 10000 -t p c
# python ~/scripts/python/topk_taxon.py -i row_sum_sanitized_summary_otu_table_mc2_w_tax_no_pynast_failures.biom.txt-rep_set_tax_assign.txt -s 0

def main(args):

    parser = argparse.ArgumentParser(description='Print top K most frequent OTUs (Taxonomies): these OTUs are most observed across samples .')
    parser.add_argument('-i', '--infile')
    parser.add_argument('-o', '--outfile', default="topk_taxon_out.txt")  # output filename
    parser.add_argument('-k', '--topk', default=5, type=int) # details for top k otus/taxonomies
    parser.add_argument('-t', '--taxons' , nargs='*')  # creates a list of items (taxons to be printed in result)
    parser.add_argument('-g', '--groupsfile')  # file containing the group where each sample belongs
    parser.add_argument('-c', '--insertcolumn', default=0, type=int) # 0-index based column where the taxonomy ids are to be printed
    parser.add_argument('-s', '--sortbycolumn', default=-1, type=int) # 0-index based column used for sorting
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    
    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile  = args.infile
    outfile = args.outfile
    topk = args.topk
    insertcolumn = args.insertcolumn
    delim = args.delimiter
    sortbycol = args.sortbycolumn

    taxon_to_print = ['p' , 'c']    
    if args.taxons != None:
        taxon_to_print = args.taxons
    
    if args.groupsfile != None:
        group_file = args.groupsfile



    # k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae; g__Kaistobacter; s__
    taxonomy_index = {'k' : 0, 'p' : 1, 'c' : 2, 'o' : 3, 'f' : 4, 'g' : 5, 's' : 6}
    # list of tuples, ordered as per the index/value of the dictionary
    sorted_taxonomy_index = sorted(taxonomy_index.iteritems(), key=operator.itemgetter(1))


    d_list=list()
    comment_list=list()

 
    for l in open(infile):
        if '#' not in l:
            d_list.append(l.strip())
        else:
            comment_list.append(l.strip())

    # descending sort as per the last column
    d_list.sort(key = lambda line: float(line.split(delim)[sortbycol]), reverse=True)

    topk_taxonom_otus_list = list()
    # add appropriate column labels to headers
    for c in comment_list:
        # walk through ordered tuple, to account for when taxon_to_print are not in correct order
        taxon_info = []
        for tup in sorted_taxonomy_index:
            if tup[0] in taxon_to_print:
                taxon_info.append(tup[0]) 
        new = '#' + delim.join(taxon_info)
        cont = c.split(delim)
        cont.insert(insertcolumn, new)
        topk_taxonom_otus_list.append(delim.join(cont))

    # print required info for topk
    for l in d_list[:topk]:
         #print l
         # split contents of line
         array = l.split(delim)
         ind=''
         taxon_info = list()
         # the taxons have already been sanitized to account for Unassigned or incomplete taxon names
         ind = [ i for i, s in enumerate(array) if re.search('; |__', s) ][0]
         # list containing taxon and its id e.g. c__abc , p__xyz
         #taxon_id_var = [l.split('__')[1] for l in array[ind].split('; ')]
         # dict containing taxon and its id e.g. c:abc , p:xyz
         tax_id_list = [i.lstrip() for i in array[ind].split(';')] # needs to be checked if it wud b same for '; '
         taxon_id_var_dict = list_to_dict(tax_id_list, delim='__')
         # walk through ordered tuple, to account for when taxon_to_print are not in correct order
         for tup in sorted_taxonomy_index:
             if tup[0] in taxon_to_print and tup[0] in taxon_id_var_dict: # check if it also exists in the taxon_id_var
                 #taxon_info.append(taxon_id_var[taxonomy_index[tup[0]]])
                 taxon_info.append(taxon_id_var_dict[tup[0]])
         array.insert(insertcolumn , delim.join(taxon_info))
         topk_taxonom_otus_list.append(delim.join(array))

    writeLIST_to_file(topk_taxonom_otus_list, outfile)

  

if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

