import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse


# 0 based index
#usage: python ~/scripts/python/row_stats.py -i sanitized_summary_otu_table_mc2_w_tax_no_pynast_failures.biom.txt-rep_set_tax_assign.txt -o row_sum_sanitized_summary_otu_table_mc2_w_tax_no_pynast_failures.biom.txt-rep_set_tax_assign.txt -s 2
#python ~/scripts/python/row_stats.py -i sanitized_summary_otu_table_mc2_w_tax_no_pynast_failures.biom.txt-rep_set_tax_assign.txt -o row_sum_sanitized_summary_otu_table_mc2_w_tax_no_pynast_failures.biom.txt-rep_set_tax_assign.txt -s 2 -c 0

def main(args):

    parser = argparse.ArgumentParser(description='(Row) sum the number of seqs. of the same OTU(taxonomy) across samples')
    parser.add_argument('-i', '--infile')
    parser.add_argument('-o', '--outfile', default="out.txt")  # output filename
    parser.add_argument('-s', '--startcolum', default=2, type=int) # start summing from start_colum,  0 based index
    parser.add_argument('-e', '--endcolum', default=99999, type=int) # end summing at (inclusive) end_colum,  0 based index
    parser.add_argument('-g', '--groupsfile')  # file containing the group where each sample belongs
    parser.add_argument('-c', '--insertcolumn', default=99999, type=int) # 0-index based column where the row_sums is to be printed
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file

    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile  = args.infile
    outpfile = args.outfile
    start_col = args.startcolum
    end_col = args.endcolum
    insertcolumn = args.insertcolumn
    delim = args.delimiter
       
    if args.groupsfile != None:
        group_file = args.groupsfile

    lines = read_file(infile)

    row_sum = '#Row_sum'
    if args.endcolum != None:
        end_col = args.endcolum + 1 # the [start:end] needs the the end value to be exclusive:1 more than last column to be summed
    
    outfile = open(outpfile, 'w')
    for l in lines:
        cont = l.strip().split(delim)
        if '#' in l:
            pass
        else:
            cont_se = list()
            cont_se = [ float(x) for x in cont[start_col:end_col] ]
            row_sum = int(sum(cont_se))
        cont.insert(insertcolumn, str(row_sum))
        outfile.write ("%s\n" % delim.join(cont))
    outfile.close()


if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

